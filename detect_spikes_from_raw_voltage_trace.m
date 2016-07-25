% Take voltage traces (saved as AD counts in .mat files as outputs of
% extract_Nlg_data.m), detect and extract potential spikes, and save as
% Neuralynx .ntt files for spike sorting in SpikeSort3D.
% 7/9/2016, Wujie Zhang
% Last updated, 7/25/2016, Wujie Zhang

%%
% Input and ouput paths, options, and paramters
voltage_trace_data_folders={'D:\Wujie\Analysis\19 Debugging and testing Nlg2Nlx code, 160629\test exports and imports\'}; % each cell is a folder where voltage traces are saved (as AD counts in .mat files as outputs of extract_Nlg_data.m)
output_folders={'D:\Wujie\Analysis\19 Debugging and testing Nlg2Nlx code, 160629\test exports and imports\'}; % each cell is a folder where the outputs of the code (time stamps and waveforms of potential spikes, in the Nlx .ntt format) will be saved, corresponding to one of the folders in voltage_trace_data_folders
output_spike_file_name_prefixes={''}; % the output file names will be a prefix followed by "TT#.ntt" where "#" is the tetrode number; can leave as an empty string; each cell here contains the prefix for one of the folders in output_folders

save_options_and_parameters=1; % 0: don't save options and paramters in a .mat file ; 1: save

filter_cutoff_frequencies=[600 6000]; % the lower and upper cut-off frequencies in Hz; the raw voltage traces will be band-pass filtered to retain only the frequency components relevant for spike detection and sorting

% Potential spikes are detected as the filtered voltage trace crosses above a threshold
hard_or_adaptive_spike_threshold=2; % 1: manually set a threshold for all channels; 2: adaptively set a threshold as a multiple of the estimated standard deviation of the noise in the voltage trace
hard_spike_threshold=[60 60 60 60]; % the threshold(s) in uV, if hard_or_adaptive_spike_threshold is 1; can either enter one number, which will be used as the threshold for all channels, or a vector with one threshold for each electrode bundle
adaptive_spike_threshold_factor=3; % the threshold is this number times the estimated noise standard deviation, if hard_or_adaptive_spike_threshold is 2; Rey et al. (2015, Brain Res Bull) recommends 3 to 5
plot_voltage_traces_with_threshold=0; % whether or not to plot some of the filtered voltage traces with the spike threshold; note: this will pause the processing and require user input to continue

spike_extraction_window_length=[-7 24]; % the first/second element indicates how many samples before/after the spike peak to extract as the spike waveform; currently (7/20/2016) the Neuralynx .ntt file needs [-7 24] for a total of 32 samples
min_separation_between_spike_peaks=abs(spike_extraction_window_length(1)); % the sample indices of the peaks of two consecutive spikes must be larger than this number; the peaks can be detected on the same or different channels

num_channels=16; % total number of recording channels, including inactive ones
channels_per_electrode_bundle=4; % the number of channels per electrode bundle, eg. 4 if using tetrodes
%%
for voltage_trace_data_folder_i=1:length(voltage_trace_data_folders) % for each of the voltage trace folders
    voltage_trace_data_folder=voltage_trace_data_folders{voltage_trace_data_folder_i};
    output_folder=output_folders{voltage_trace_data_folder_i};
    output_spike_file_name_prefix=output_spike_file_name_prefixes{voltage_trace_data_folder_i};
    disp(['Processing the data in "' voltage_trace_data_folder '"...'])
    %%
    num_electrode_bundle=num_channels/channels_per_electrode_bundle; % number of electrode bundles, eg. tetrodes
    if plot_voltage_traces_with_threshold
        fig_threshold_check=figure;
    end
    voltage_data_list=dir([voltage_trace_data_folder 'CSC*.mat']); % find all the .mat format voltage data files in the folder
    active_channels=zeros(1,length(voltage_data_list)); % each element will be a number between 1 and the total number of channels, indicating the number of each active channel
    for file_i=1:length(voltage_data_list)
        active_channels(file_i)=str2double(voltage_data_list(file_i).name(4:end-4))+1; % here we number the channels starting from 1, unlike the file names which follow the Neuralynx convention of starting from 0
    end
    if channels_per_electrode_bundle==4
        bundle_name='Tetrode';
    else
        bundle_name='Electrode bundle';
    end
    
    first_file_loaded=1;
    if ~exist(output_folder,'dir'); % make the output folder if it doesn't already exist
        mkdir(output_folder);
    end
    for electrode_bundle_i=1:num_electrode_bundle % for each of the electrode bundles, eg. tetrodes
        active_channels_on_current_bundle=intersect(active_channels,(electrode_bundle_i-1)*channels_per_electrode_bundle+1:electrode_bundle_i*channels_per_electrode_bundle); % find the active channels on the current electrode bundle
        if isempty(active_channels_on_current_bundle)
            continue
        end
        num_channels_on_current_bundle=length(active_channels_on_current_bundle);
        
        % Filtering the voltage traces
        disp(['Filtering voltage traces from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        for current_bundle_channel_i=1:num_channels_on_current_bundle % for each of the active channels on this electrode bundle
            CSC_num=active_channels_on_current_bundle(current_bundle_channel_i)-1;
            load(fullfile(voltage_trace_data_folder,['CSC' num2str(CSC_num) '.mat']))
            if first_file_loaded==1 % this is true for the first file loaded
                first_file_loaded=0; % this is set to zero after the first file has already been loaded
                num_samples_per_channel=length(AD_count_int16);
                start_recording_during_session=1+find(diff(indices_of_first_samples)~=samples_per_channel_per_file); % the indices of the Nlg .DAT files where recording was started after being stopped mid-session; note these indices start from 1, unlike the actual file names of the Nlg .DAT files which start from 000; this line finds them by finding the files where recording was stopped, which have less than the normal number of samples as a result of processing by extract_Nlg_data.m
                discontinuous_start_sample_indices=indices_of_first_samples(start_recording_during_session); % for each of these samples, there was a large time gap between it and the immediately previous sample
                recording_start_end_sample_indices=[1 discontinuous_start_sample_indices num_samples_per_channel+1]; % filtering will break the voltage traces into chunks that go from recording_start_end_sample_indices(i) to recording_start_end_sample_indices(i+1)-1
            end
            if current_bundle_channel_i==1
                filtered_voltage_traces=zeros(num_channels_on_current_bundle,num_samples_per_channel); % preallocate the variable
            end
            sampling_freq=1/(sampling_period_usec/1e6); % in Hz
            voltage_trace=double(AD_count_int16)*AD_count_to_uV_factor; % convert to voltages in uV in "double" format
            [b,a]=butter(6,filter_cutoff_frequencies/(sampling_freq/2),'bandpass'); % a Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency, as required by MATLAB)
            for continuous_recording_period_i=1:length(recording_start_end_sample_indices)-1 % filtering the chunks of the voltage trace separately, because they are not continuous in time
                period_start=recording_start_end_sample_indices(continuous_recording_period_i);
                period_end=recording_start_end_sample_indices(continuous_recording_period_i+1)-1;
                filtered_voltage_traces(current_bundle_channel_i,period_start:period_end)=filtfilt(b,a,voltage_trace(period_start:period_end)); % band-pass filter the voltage traces
            end
        end
        clear AD_count_int16 voltage_trace
        
        % Detect spikes as threshold-crossing by the filtered voltage traces
        if hard_or_adaptive_spike_threshold==1 % using the manually-set hard threshold(s)
            if length(hard_spike_threshold)==1
                spike_thresholds=hard_spike_threshold*ones(num_channels_on_current_bundle,1);
            elseif length(hard_spike_threshold)==num_electrode_bundle
                spike_thresholds=hard_spike_threshold(electrode_bundle_i)*ones(num_channels_on_current_bundle,1);
            end
        elseif hard_or_adaptive_spike_threshold==2 % automatically calculate a threshold (from Quian Quiroga et al., 2004, Neural Computation)
            estimate_of_voltage_noise_std=median(abs(filtered_voltage_traces),2)/icdf('Normal',0.75,0,1);
            % Assume the fluctuations of the filtered voltage during
            % non-spiking periods (ie. the noise) is normally distributed
            % around zero. The median of the absolute values of the noise
            % (approximated here by the median of the absolute values of the
            % entire voltage trace, because spikes constitute a small
            % proportion of the voltage trace) is the 75th percentile of the
            % noise distribution. Then, this median divided by the 75th
            % percentile of the standard normal distribution (whose standard
            % deviation is 1) is the standard deviation of the noise
            % distribution.
            spike_thresholds=adaptive_spike_threshold_factor*estimate_of_voltage_noise_std;
        end
        
        if plot_voltage_traces_with_threshold
            figure(fig_threshold_check)
            num_samples_to_plot=500; % number of samples to plot for each plotting window
            for window_i=1:50 % plot some windows sequentially
                clf
                samples_to_plot=(window_i-1)*num_samples_to_plot+1:window_i*num_samples_to_plot;
                times_to_plot=(1:length(samples_to_plot))*(sampling_period_usec/1000); % in ms
                for plot_channel_i=1:num_channels_on_current_bundle
                    subplot(num_channels_on_current_bundle,1,plot_channel_i)
                    hold on
                    plot(times_to_plot,filtered_voltage_traces(plot_channel_i,samples_to_plot),'k')
                    plot(times_to_plot,repmat(spike_thresholds(plot_channel_i),1,length(times_to_plot)),'r--')
                    ylabel(['Channel ' num2str(active_channels_on_current_bundle(plot_channel_i)) ' voltage (uV)'])
                    xlim(times_to_plot([1 end]))
                    if plot_channel_i==1
                        title([bundle_name ' ' num2str(electrode_bundle_i)])
                    end
                end
                xlabel('Time (ms)')
                legend('Filtered voltage trace','Spike threshold')
                stop_plotting=input('Enter anything to stop plotting: ','s');
                if any(stop_plotting)
                    break
                end
            end
        end
        
        disp(['Detecting spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        voltage_peaks=zeros(1,num_samples_per_channel); % a vector that will store, for each time point, the highest above-threshold peak across all channels, where each peak is detected on a single channel and its height is normalized to the channel's threshold
        for current_bundle_channel_i=1:num_channels_on_current_bundle % for each of the active channels on the current electrode bundle
            [peak_heights,sample_indices_of_peaks]=findpeaks(filtered_voltage_traces(current_bundle_channel_i,:),'MinPeakHeight',spike_thresholds(current_bundle_channel_i)); % find all the peaks in the filtered voltage trace that are also above threshold
            peak_heights=peak_heights/spike_thresholds(current_bundle_channel_i); % normalize the peak heights by the threshold for that channel
            higher_peaks=voltage_peaks(sample_indices_of_peaks)<peak_heights; % compare the heights of the peaks already stored with the ones detected on the current channel, and use logical indexing to express where the current channel has the higher peaks
            voltage_peaks(sample_indices_of_peaks(higher_peaks))=peak_heights(higher_peaks); % replace the stored peak heights with the ones detected on the current channel, wherever the current channel has higher peaks
        end
        [~,sample_indices_of_peaks]=findpeaks(voltage_peaks,'MinPeakDistance',min_separation_between_spike_peaks); % find peaks that are separated by a minimum distance; if two peaks are within that distance, only find the higher peak
        sample_indices_of_peaks(sample_indices_of_peaks+spike_extraction_window_length(1)<1 | sample_indices_of_peaks+spike_extraction_window_length(2)>num_samples_per_channel)=[]; % delete the peaks that would lead to spike waveform windows that extend beyond the recording time
        clear voltage_peaks
        
        % Extracting the spike waveforms
        disp(['Extracting spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        num_spikes=length(sample_indices_of_peaks);
        spike_waveforms=zeros(32,4,num_spikes);
        timestamps_usec=round(get_timestamps_for_Nlg_voltage_samples(sample_indices_of_peaks,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)); % the time stamps of all the spike peaks, rounded to integer microseconds
        for spike_i=1:num_spikes
            channel_indices_for_current_bundle=active_channels_on_current_bundle-(electrode_bundle_i-1)*channels_per_electrode_bundle; % eg. converting channels 5, 6, 7, and 8 (the four electrodes on tetrode 2) into indices 1, 2, 3, and 4
            spike_waveforms(:,channel_indices_for_current_bundle,spike_i)=round(filtered_voltage_traces(:,sample_indices_of_peaks(spike_i)+spike_extraction_window_length(1):sample_indices_of_peaks(spike_i)+spike_extraction_window_length(2))).'; % save the waveforms of the current spike from all active channels, rounding to integer uV
        end
        
        % Saving the time stamps and waveforms of spikes in Nlx .ntt format
        disp(['Saving spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        % The following works with the current version of Mat2NlxSpike (Version
        % 6.0.0) as of this writing (7/13/2016).
        NTT_file_name=[output_spike_file_name_prefix 'TT' num2str(electrode_bundle_i) '.ntt'];
        file_name_to_save=fullfile(output_folder,NTT_file_name);
        AppendToFileFlag=0; % save new file, or overwrite existing file, but do not append to existing file
        ExportMode=1; % export all spikes
        ExportModeVector=[]; % don't need this when exporting all spikes
        FieldSelectionFlags=[1 0 1 0 1 1]; % whether or not to export: time stamps, spike channel numbers, the clusters that spikes are assigned to, spike features, spike waveforms, and header
        cluster_numbers=zeros(1,length(timestamps_usec)); % zero means that the corresponding spike is not assigned to any cluster; we have to do this because Mat2NlxSpike 6.0.0 has a bug: if this input is omitted, all spikes are incorrectly assigned to clusters
        NTT_header={'######## Neuralynx Data File Header';'-ADMaxValue 500';'-ADBitVolts 1.0 1.0 1.0 1.0'}; % "-ADMaxValue 500" ensures that when "shift + right click" to preview waveforms in SpikeSort3D, the voltage axis has reasonable scales
        Mat2NlxSpike(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,timestamps_usec,cluster_numbers,spike_waveforms,NTT_header)
        % Note that if the time stamp and waveform inputs to Mat2NlxSpike
        % (6.0.0) are in certain integer formats, Mat2NlxSpike cannot correctly
        % save them, even though the .NTT file saves the waveforms and time
        % stamps in integer data formats. Thus, "timestamps_usec" and
        % "spike_waveforms" here are in the "double" format.
    end
    disp(['Finished saving spikes from all voltage traces in "' voltage_trace_data_folder '".'])
    %%
    % Save the options and parameters that were used
    if save_options_and_parameters
        file_name_to_save=fullfile(output_folder,['detect_spikes_from_raw_voltage_trace_paramters_' date '.mat']);
        date_time_of_processing=datetime; % the date and time when this code was run
        save(file_name_to_save,'voltage_trace_data_folder','output_folder','output_spike_file_name_prefix','filter_cutoff_frequencies','hard_or_adaptive_spike_threshold','hard_spike_threshold','adaptive_spike_threshold_factor','spike_extraction_window_length','min_separation_between_spike_peaks','num_channels','channels_per_electrode_bundle','spike_thresholds','date_time_of_processing')
    end
end
