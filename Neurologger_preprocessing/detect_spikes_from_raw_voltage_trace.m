% Take voltage traces (saved as AD counts in .mat files as outputs of
% extract_Nlg_data.m), detect and extract potential spikes, and save as
% Neuralynx .ntt files for spike sorting in SpikeSort3D.
% 7/9/2016, Wujie Zhang
last_code_update='5/16/2017, Wujie Zhang'; % identifies the version of the code
%%
% Input and ouput paths, options, and paramters
voltage_trace_data_folders={'H:\Wujie\Data\two bat recording 1\20170304\communication\bat59813_YL\neural\' 'H:\Wujie\Data\two bat recording 1\20170304\communication\bat60141_TC\neural\'}; % each cell is a folder where voltage traces are saved (as AD counts in .mat files as outputs of extract_Nlg_data.m)
output_folders=voltage_trace_data_folders; % each cell is a folder where the outputs of the code (time stamps and waveforms of potential spikes, in the Nlx .ntt format) will be saved, corresponding to one of the folders in voltage_trace_data_folders

save_options_and_parameters=1; % 0: don't save options and paramters in a .mat file ; 1: save
output_spike_file_name_prefix=''; % the output file names will be a prefix followed by "TT#.ntt" where "#" is the tetrode number; can leave as an empty string
filter_cutoff_frequencies=[600 6000]; % the lower and upper cut-off frequencies in Hz; the raw voltage traces will be band-pass filtered to retain only the frequency components relevant for spike detection and sorting

% Potential spikes are detected as the filtered voltage trace crosses above a threshold
manual_or_automatic_spike_threshold=2; % 1: manually set a threshold for all channels; 2: automatically set a threshold as a multiple of the estimated standard deviation of the noise in the voltage trace
manual_spike_threshold=[40 40 40 40]; % the threshold(s) in uV, if manual_or_automatic_spike_threshold is 1; can either enter one number, which will be used as the threshold for all channels, or a vector with one threshold for each electrode bundle
automatic_spike_threshold_factor=3; % the threshold is this number times the estimated noise standard deviation, if manual_or_automatic_spike_threshold is 2; Rey et al. (2015, Brain Res Bull) recommends 3 to 5
absolute_or_threshold_normalized_peak_heights=1; % when comparing the heights of voltage peaks on different channels of the same electrode bundle, whether to use 1: the absolute peak heights; or 2: peak heights normalized by the spike thresholds of the respective channels; this only makes a difference if using the automatic thresholds
plot_voltage_traces_with_threshold=0; % whether or not to plot some of the filtered voltage traces with the spike threshold; note: this will pause the processing and require user input to continue

% Comparison of the waveforms of detected spikes with a library of accepted
% waveforms; adapted from Michael Yartsev
compare_detected_waveforms_with_library_or_typical_waveforms=2; % 1: for each potential spike, calculate the correlations between its waveform and the waveforms from a library of accepted spikes, and reject the waveforms with low corelations; 2: calculate correlations with only two typical waveforms from the library (one broad and one narrow spike); 0: don't do this step
waveform_library_file='D:\Wujie\Scripts\Code package for prepocessing Nlg data\library_of_acceptable_spike_shapes_peak_aligned_at_8.mat'; % the library of acceptable waveforms compiled by Michael Yartsev
lowest_acceptable_correlation=0.5; % for each potential spike, if the highest correlation is below this value, then reject that potential spike; Michael recommends 0.95 if calculating correlations with the entire library

% Check if all electrode bundles (tetrodes) detect spikes at the same
% times; adapted from Michael Yartsev
check_spike_coincidence_across_electrode_bundles=1; % 1: if a potential spike appears at the same time across all electode bundles (tetrodes), then delete it because it's likely an artifact; 0: don't do this
maximum_coincidence_interval_usec=50; % if the maximum difference between the times of spike detection from all electode bundles (tetrodes) is less than or equal to this value, then the spike is deleted on all electode bundles (tetrodes); in microsec

spike_extraction_window_length=[-7 24]; % the first/second element indicates how many samples before/after the spike peak to extract as the spike waveform; currently (7/20/2016) the Neuralynx .ntt file needs [-7 24] for a total of 32 samples
min_separation_between_spike_peaks=abs(spike_extraction_window_length(1)); % the sample indices of the peaks of two consecutive spikes must be larger than this number; the peaks can be detected on the same or different channels

num_channels=16; % total number of recording channels, including inactive ones
channels_per_electrode_bundle=4; % the number of channels per electrode bundle, eg. 4 if using tetrodes
%%
for voltage_trace_data_folder_i=1:length(voltage_trace_data_folders) % for each of the voltage trace folders
    voltage_trace_data_folder=voltage_trace_data_folders{voltage_trace_data_folder_i};
    output_folder=output_folders{voltage_trace_data_folder_i};
    disp(['Processing the data in "' voltage_trace_data_folder '"...'])
    %%
    num_electrode_bundle=num_channels/channels_per_electrode_bundle; % number of electrode bundles, eg. tetrodes
    if plot_voltage_traces_with_threshold
        fig_threshold_check=figure;
    end
    voltage_data_list=dir(fullfile(voltage_trace_data_folder,'CSC*.mat')); % find all the .mat format voltage data files in the folder
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
    if ~exist(output_folder,'dir') % make the output folder if it doesn't already exist
        mkdir(output_folder);
    end
    spike_thresholds_all_channels=nan(channels_per_electrode_bundle*num_electrode_bundle,1);
    timestamps_usec_all_electrode_bundles=cell(num_electrode_bundle,1);
    spike_waveforms_all_electrode_bundles=cell(num_electrode_bundle,1);
    num_spikes_all_electrode_bundles=nan(num_electrode_bundle,1);
    
    for electrode_bundle_i=1:num_electrode_bundle % for each of the electrode bundles, eg. tetrodes
        active_channels_on_current_bundle=intersect(active_channels,(electrode_bundle_i-1)*channels_per_electrode_bundle+1:electrode_bundle_i*channels_per_electrode_bundle); % find the active channels on the current electrode bundle
        if isempty(active_channels_on_current_bundle)
            continue
        end
        num_channels_on_current_bundle=length(active_channels_on_current_bundle);
        
        %%
        % Filtering the voltage traces
        disp(['Filtering voltage traces from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        for current_bundle_channel_i=1:num_channels_on_current_bundle % for each of the active channels on this electrode bundle
            CSC_num=active_channels_on_current_bundle(current_bundle_channel_i)-1;
            load(fullfile(voltage_trace_data_folder,['CSC' num2str(CSC_num) '.mat']))
            if first_file_loaded==1 % this is true for the first file loaded
                first_file_loaded=0; % this is set to zero after the first file has already been loaded
                num_samples_per_channel=length(AD_count_int16);
                start_recording_during_session=1+find(diff(indices_of_first_samples)~=samples_per_channel_per_file); % the indices of the Nlg .DAT files where recording was started after being stopped mid-session; note these indices start from 1, unlike the actual file names of the Nlg .DAT files which start from 000; this line finds them by finding the files where recording was stopped, which have less than the normal number of samples, because extract_Nlg_data.m deleted the samples that do not contain actual recorded data
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
        
        %%
        % Detect spikes as threshold-crossing by the filtered voltage traces
        if manual_or_automatic_spike_threshold==1 % using the manually-set threshold(s)
            if length(manual_spike_threshold)==1
                spike_thresholds=manual_spike_threshold*ones(num_channels_on_current_bundle,1);
            elseif length(manual_spike_threshold)==num_electrode_bundle
                spike_thresholds=manual_spike_threshold(electrode_bundle_i)*ones(num_channels_on_current_bundle,1);
            end
        elseif manual_or_automatic_spike_threshold==2 % automatically calculate a threshold (modified from Quian Quiroga et al., 2004, Neural Computation)
            estimate_of_voltage_noise_std=(quantile(filtered_voltage_traces,0.75,2)-median(filtered_voltage_traces,2))/icdf('Normal',0.75,0,1);
            % Assume the fluctuations of the filtered voltage during
            % non-spiking periods (ie. the noise) is normally distributed.
            % We approximate the median and the 75th percentile of the
            % noise voltage by the median and 75th percentile of the entire
            % voltage trace, because spikes constitute a small proportion
            % of the voltage trace. Then the difference between the 75th
            % percentile and the median, divided by the 75th percentile of
            % the standard normal distribution (whose mean is 0 and
            % standard deviation is 1) is the standard deviation of the
            % noise distribution. Because the Neurologger voltage
            % measurements may sometimes have an offset (eg. all voltages
            % are shifted slightly up from their true values), here we
            % don't assume the noise distribution is centered around zero,
            % unlike Quian Quiroga et al. (2004).
            spike_thresholds=automatic_spike_threshold_factor*estimate_of_voltage_noise_std;
        end
        spike_thresholds_all_channels(active_channels_on_current_bundle)=spike_thresholds;
        
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
            if absolute_or_threshold_normalized_peak_heights==2
                peak_heights=peak_heights/spike_thresholds(current_bundle_channel_i); % normalize the peak heights by the threshold for that channel
            end
            higher_peaks=voltage_peaks(sample_indices_of_peaks)<peak_heights; % compare the heights of the peaks already stored with the ones detected on the current channel, and use logical indexing to express where the current channel has the higher peaks
            voltage_peaks(sample_indices_of_peaks(higher_peaks))=peak_heights(higher_peaks); % replace the stored peak heights with the ones detected on the current channel, wherever the current channel has higher peaks
        end
        [~,sample_indices_of_peaks]=findpeaks(voltage_peaks,'MinPeakDistance',min_separation_between_spike_peaks); % find peaks that are separated by a minimum distance; if two peaks are within that distance, only find the higher peak
        sample_indices_of_peaks(sample_indices_of_peaks+spike_extraction_window_length(1)<1 | sample_indices_of_peaks+spike_extraction_window_length(2)>num_samples_per_channel)=[]; % delete the peaks that would lead to spike waveform windows that extend beyond the recording time
        clear voltage_peaks
        
        %%
        % Extracting the spike waveforms
        disp(['Extracting spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
        num_spikes=length(sample_indices_of_peaks);
        spike_waveforms=zeros(32,4,num_spikes);
        timestamps_usec=round(get_timestamps_for_Nlg_voltage_samples(sample_indices_of_peaks,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)); % the time stamps of all the spike peaks, rounded to integer microseconds; note that these are the time stamps of the last channel on this electrode bundle, which differ from the time stamps on the other channels of this electrode bundle by a few sampling periods of the Nlg AD converter
        for spike_i=1:num_spikes
            channel_indices_for_current_bundle=active_channels_on_current_bundle-(electrode_bundle_i-1)*channels_per_electrode_bundle; % eg. converting channels 5, 6, 7, and 8 (the four electrodes on tetrode 2) into indices 1, 2, 3, and 4
            spike_waveforms(:,channel_indices_for_current_bundle,spike_i)=filtered_voltage_traces(:,sample_indices_of_peaks(spike_i)+spike_extraction_window_length(1):sample_indices_of_peaks(spike_i)+spike_extraction_window_length(2)).'; % save the waveforms of the current spike from all active channels (units are uV)
        end
        
        %%
        % Comparison of the waveforms of detected spikes with a library of
        % accepted waveforms; adapted from Michael Yartsev
        if compare_detected_waveforms_with_library_or_typical_waveforms
            load(waveform_library_file) % this contains the variable "library_of_acceptable_spike_shapes" which is 183 rows (the different waveforms) by 32 columns (the 32 voltage samples of each waveforms); the waveforms are all peak-normalized
            logical_indices_of_accepted_spikes=false(size(spike_waveforms,3),1); % a given element of this variable will be 1 (or 0) if the corresponding spike will be accepted (or rejected)
            disp(['Comparing detected spike waveforms from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) ' with waveforms from library...'])
            for spike_i=1:size(spike_waveforms,3)
                [~,channel_with_highest_peak]=max(max(spike_waveforms(:,:,spike_i),[],1)); % find the channel with the highest peak
                waveform_to_test=spike_waveforms(:,channel_with_highest_peak,spike_i);
                if compare_detected_waveforms_with_library_or_typical_waveforms==1
                    max_correlation=max([corr(waveform_to_test(2:end-1),library_of_acceptable_spike_shapes(:,2:end-1).') corr(waveform_to_test(2:end-1),library_of_acceptable_spike_shapes(:,1:end-2).') corr(waveform_to_test(2:end-1),library_of_acceptable_spike_shapes(:,3:end).')]); % for each detected spike, Michael calculates correlations after shifting the waveforms in the library by -1, 0, and 1 samples, and look for the maximum correlation
                elseif compare_detected_waveforms_with_library_or_typical_waveforms==2
                    max_correlation=max(corr(waveform_to_test,library_of_acceptable_spike_shapes([15 176],:).')); % the 15th and 176th waveforms in the library are typical broad and narrow spikes, respectively; almost all other waveforms in the library have at least 0.95 correlation with at least one of these two waveforms
                end
                logical_indices_of_accepted_spikes(spike_i)=max_correlation>=lowest_acceptable_correlation;
            end
            if sum(logical_indices_of_accepted_spikes)<length(timestamps_usec)
                disp(['Rejected ' num2str(length(timestamps_usec)-sum(logical_indices_of_accepted_spikes)) ' potential spikes after comparison with waveform library...'])
            end
            spike_waveforms=spike_waveforms(:,:,logical_indices_of_accepted_spikes);
            timestamps_usec=timestamps_usec(logical_indices_of_accepted_spikes);
        end
        
        %%
        timestamps_usec_all_electrode_bundles{electrode_bundle_i}=timestamps_usec;
        spike_waveforms_all_electrode_bundles{electrode_bundle_i}=spike_waveforms;
        num_spikes_all_electrode_bundles(electrode_bundle_i)=length(timestamps_usec);
    end
    
    %%
    % If all electrode bundles (tetrodes) detect a spike at the same time
    % (ie. within a coincidence time window), it's possible that that's an
    % artifact; here we check for and reject all such spikes; adapted from
    % Michael Yartsev
    % Note: currently, this code does not rigorously deal with uncommon
    % cases when more than 1 spike occurs within the coincidence time
    % window on the SAME electrode bundle
    if check_spike_coincidence_across_electrode_bundles
        disp(['Checking for spikes detected at the same time on all ' lower(bundle_name) 's...'])
        indices_active_bundles=find(~isnan(num_spikes_all_electrode_bundles)).'; % indices of the active electrode bundles (tetrodes)
        if length(indices_active_bundles)>1
            [~,index_bundle_with_minimum_num_spikes]=min(num_spikes_all_electrode_bundles); % index of the electrode bundle with the smallest number of spikes; the number of spikes on inactive bundles is NaN, which does not count as smaller than actual numbers here
            spike_indices_to_delete=zeros(num_electrode_bundle,length(timestamps_usec_all_electrode_bundles{index_bundle_with_minimum_num_spikes})); % each row is an electrode bundle, the number of columns is the maximum number of spikes that can be detected at the same time on all electrode bundles
            num_detected_coincidences=0;
            for bundle_spike_i=1:length(timestamps_usec_all_electrode_bundles{index_bundle_with_minimum_num_spikes}) % for each of the spikes from the electrode bundle with the smallest number of spikes
                go_to_next_spike=0;
                current_spike_times_to_consider_for_coincidence=nan(num_electrode_bundle,1); % for the current spike that is possibly occuring at the same time on all electrode bundles, this variable will be filled with the spike times on the different electrode bundles
                current_spike_times_to_consider_for_coincidence(index_bundle_with_minimum_num_spikes)=timestamps_usec_all_electrode_bundles{index_bundle_with_minimum_num_spikes}(bundle_spike_i);
                current_spike_indices_to_consider_for_coincidence=nan(num_electrode_bundle,1); % the indices within each electrode bundle, of the spikes that are currently being considered as occuring at the same time
                current_spike_indices_to_consider_for_coincidence(index_bundle_with_minimum_num_spikes)=bundle_spike_i;
                for electrode_bundle_i=indices_active_bundles(~ismember(indices_active_bundles,index_bundle_with_minimum_num_spikes)) % for each of the other electrode bundles
                    intervals_between_spikes_on_two_bundles=abs(timestamps_usec_all_electrode_bundles{electrode_bundle_i}-current_spike_times_to_consider_for_coincidence(index_bundle_with_minimum_num_spikes)); % the intervals between all the spikes on electrode_bundle_i and the current spike on the electrode bundle with the smallest number of spikes
                    [min_interval,index_closest_spike]=min(intervals_between_spikes_on_two_bundles); % the spike on electrode_bundle_i that is closest
                    if min_interval<=maximum_coincidence_interval_usec
                        current_spike_times_to_consider_for_coincidence(electrode_bundle_i)=timestamps_usec_all_electrode_bundles{electrode_bundle_i}(index_closest_spike);
                        current_spike_indices_to_consider_for_coincidence(electrode_bundle_i)=index_closest_spike;
                    else
                        go_to_next_spike=1;
                        break
                    end
                end
                if go_to_next_spike==0 && range(current_spike_times_to_consider_for_coincidence)<=maximum_coincidence_interval_usec % check if the spike times on the other electode bundles are all within the coincidence time window; the "range" function disregards NaNs, which would be on inactive electrode bundles
                    num_detected_coincidences=num_detected_coincidences+1;
                    spike_indices_to_delete(:,num_detected_coincidences)=current_spike_indices_to_consider_for_coincidence;
                end
            end
            for electrode_bundle_i=indices_active_bundles
                timestamps_usec_all_electrode_bundles{electrode_bundle_i}(spike_indices_to_delete(electrode_bundle_i,1:num_detected_coincidences))=[]; % delete the rejected spikes
                spike_waveforms_all_electrode_bundles{electrode_bundle_i}(:,:,spike_indices_to_delete(electrode_bundle_i,1:num_detected_coincidences))=[];
            end
            disp(['Deleted ' num2str(num_detected_coincidences) ' potential spikes that were detected at the same time on all ' lower(bundle_name) 's...'])
        else
            disp(['There are less than two ' lower(bundle_name) 's; skipping spike coincidence detection across ' lower(bundle_name) 's...'])
        end
    end
    
    %%
    for electrode_bundle_i=1:num_electrode_bundle % for each of the electrode bundles, eg. tetrodes
        if ~isempty(timestamps_usec_all_electrode_bundles{electrode_bundle_i}) % check if the electrode bundle has no spikes (eg. if all its channels are inactive)
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
            cluster_numbers=zeros(1,length(timestamps_usec_all_electrode_bundles{electrode_bundle_i})); % zero means that the corresponding spike is not assigned to any cluster; we have to do this because Mat2NlxSpike 6.0.0 has a bug: if this input is omitted, all spikes are incorrectly assigned to clusters
            
            % To pan the 2D view in SpikeSort3D, "-ADMaxValue 32767" must
            % be in the header; to ensure that the voltage axis has
            % reasonable scales when using "shift + right click" to preview
            % waveforms in SpikeSort3D, we have to scale the waveforms up
            % as follows
            ADBitVolts=500/(32767*10^6); % 32767 is the largest number that can be represented as a signed 16-bit integer, which SpikeSort3D requires to be mapped to the upper bound of the voltage data, which we have picked to be 500 uV here, so that the voltage axis during preview has reasonable scales
            waveforms_to_export_AD_counts=round(spike_waveforms_all_electrode_bundles{electrode_bundle_i}/(ADBitVolts*10^6)); % the AD counts here, multiplied by the "ADBitVolts" factor above, equals voltages in V
            NTT_header={'######## Neuralynx Data File Header';'-ADMaxValue 32767';['-ADBitVolts ' sprintf('%.24f',ADBitVolts) ' ' sprintf('%.24f',ADBitVolts) ' ' sprintf('%.24f',ADBitVolts) ' ' sprintf('%.24f',ADBitVolts)]}; % the string repetitions of "ADBitVolts" are for the different channels; the number is saved into the header with 24 digits after the decimal point, as the Neuralynx recording system does; note that this does no round to the 24th digit after the decimal point, which is OK when "ADBitVolts" is 500/(32767*10^6), but needs to be checked if a different voltage upper bound is used
            
            Mat2NlxSpike(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,timestamps_usec_all_electrode_bundles{electrode_bundle_i},cluster_numbers,waveforms_to_export_AD_counts,NTT_header)
            % Note that if the time stamp and waveform inputs to Mat2NlxSpike
            % (6.0.0) are in certain integer formats, Mat2NlxSpike cannot correctly
            % save them, even though the .NTT file saves the waveforms and time
            % stamps in integer data formats. Thus, "timestamps_usec" and
            % "spike_waveforms" here are integers in the "double" format.
        end
    end
    disp(['Finished saving spikes from all voltage traces in "' voltage_trace_data_folder '".'])
    %%
    % Save the options and parameters that were used
    if save_options_and_parameters
        file_name_to_save=fullfile(output_folder,['detect_spikes_from_raw_voltage_trace_paramters_' date '.mat']);
        date_time_of_processing=datetime; % the date and time when this code was run
        clear variables_to_save
        variables_to_save.voltage_trace_data_folder=voltage_trace_data_folder;
        variables_to_save.output_folder=output_folder;
        variables_to_save.output_spike_file_name_prefix=output_spike_file_name_prefix;
        variables_to_save.filter_cutoff_frequencies=filter_cutoff_frequencies;
        variables_to_save.manual_or_automatic_spike_threshold=manual_or_automatic_spike_threshold;
        variables_to_save.manual_spike_threshold=manual_spike_threshold;
        variables_to_save.automatic_spike_threshold_factor=automatic_spike_threshold_factor;
        variables_to_save.spike_extraction_window_length=spike_extraction_window_length;
        variables_to_save.min_separation_between_spike_peaks=min_separation_between_spike_peaks;
        variables_to_save.num_channels=num_channels;
        variables_to_save.channels_per_electrode_bundle=channels_per_electrode_bundle;
        variables_to_save.spike_thresholds_all_channels=spike_thresholds_all_channels;
        variables_to_save.absolute_or_threshold_normalized_peak_heights=absolute_or_threshold_normalized_peak_heights;
        variables_to_save.compare_detected_waveforms_with_library=compare_detected_waveforms_with_library_or_typical_waveforms;
        variables_to_save.waveform_library_file=waveform_library_file;
        variables_to_save.lowest_acceptable_correlation=lowest_acceptable_correlation;
        variables_to_save.check_spike_coincidence_across_electrode_bundles=check_spike_coincidence_across_electrode_bundles;
        variables_to_save.maximum_coincidence_interval_usec=maximum_coincidence_interval_usec;
        variables_to_save.date_time_of_processing=date_time_of_processing;
        variables_to_save.last_code_update=last_code_update;
        save(file_name_to_save,'-struct','variables_to_save')
    end
end
