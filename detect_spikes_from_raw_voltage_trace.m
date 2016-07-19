% Take voltage traces (saved in .mat files as outputs of
% extract_Nlg_data.m), detect and extract potential spikes, and save as
% Neuralynx .ntt files for spike sorting in SpikeSort3D.
% 7/9/2016, Wujie Zhang
% Beta testing version. Do not use this to seriously process data.
%%
voltage_trace_data_folder='C:\Users\phyllo\Documents\Maimon\ephys\yr2016_bat71319_robin_Nlg\neurologger_recording20160718\nlxformat\';
output_folder='C:\Users\phyllo\Documents\Maimon\ephys\Data_processed\yr2016_bat71319_robin_Nlg\20160718\';
output_spike_file_name_prefix='20160718';

filter_cutoff_frequencies=[600 6000];

% Potential spikes are detected as a filtered voltage trace crosses above a threshold
hard_or_adaptive_spike_threshold=2; % 1: manually set a threshold for all channels; 2: adaptively set a threshold as a multiple of the estimated standard deviation of the noise in the voltage trace
hard_spike_threshold=60; % use this threshold (in uV), if hard_or_adaptive_spike_threshold is 1
adaptive_spike_threshold_factor=3; % the threshold is this number times the estimated noise standard deviation, if hard_or_adaptive_spike_threshold is 2; Rey et al. (2015, Brain Res Bull) recommends 3 to 5.
plot_voltage_traces_with_threshold=0; % whether or not to plot the filtered voltage traces with the spike threshold to check the threshold

spike_extraction_window_length=[-7 24]; % the first/second element indicates how many samples before/after the spike peak to extract as the spike waveform; currently the Neuralynx .ntt file needs [-7 24] for a total of 32 samples
min_separation_between_spike_peaks=abs(spike_extraction_window_length(1)); % the sample indices of the peaks of two consecutive spikes must be larger than this number; the peaks can be detected on the same or different channels

num_channels=16;
channels_per_electrode_bundle=4; % the number of channels per electrode bundle, eg. 4 if using tetrodes
AD_count_to_uV_factor=3.3; % voltage in microvolt = AD count * AD_count_to_uV_factor

%%
num_electrode_bundle=num_channels/channels_per_electrode_bundle;
if plot_voltage_traces_with_threshold
    fig_threshold_check=figure;
end
voltage_data_list=dir([voltage_trace_data_folder 'CSC*.mat']);
active_channels=zeros(1,length(voltage_data_list));
for file_i=1:length(voltage_data_list)
    active_channels(file_i)=str2double(voltage_data_list(file_i).name(4:end-4))+1; % here we number the channels starting from 1, unlike the file names which follow the Neuralynx convention of starting from 0
end
if channels_per_electrode_bundle==4
    bundle_name='Tetrode';
else
    bundle_name='Electrode bundle';
end

first_file_loaded=1;
if ~exist(output_folder,'dir');
    mkdir(output_folder);
end
for electrode_bundle_i=1:num_electrode_bundle
    active_channels_on_current_bundle=intersect(active_channels,(electrode_bundle_i-1)*channels_per_electrode_bundle+1:electrode_bundle_i*channels_per_electrode_bundle);
    if isempty(active_channels_on_current_bundle)
        continue
    end
    
    num_channels_on_current_bundle=length(active_channels_on_current_bundle);
    disp(['Filtering voltage traces from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
    for current_bundle_channel_i=1:num_channels_on_current_bundle
        CSC_num=active_channels_on_current_bundle(current_bundle_channel_i)-1;
        load(fullfile(voltage_trace_data_folder,['CSC' num2str(CSC_num) '.mat']))
        if first_file_loaded==1
            first_file_loaded=0;
            num_samples_per_channel=length(AD_count_int16);
            start_recording_during_session=1+find(diff(indices_of_first_samples)~=samples_per_channel_per_file); % the indices of the Nlg .DAT files where recording was started after being stopped mid-session; note these indices start from 1, unlike the actual file names of the Nlg .DAT files which start from 000
            discontinuous_start_sample_indices=indices_of_first_samples(start_recording_during_session); % for each of these samples, there was a large time gap between it and the immediately previous sample
            recording_start_end_sample_indices=[1 discontinuous_start_sample_indices num_samples_per_channel+1];
        end
        if current_bundle_channel_i==1
            filtered_voltage_traces=zeros(num_channels_on_current_bundle,num_samples_per_channel);
        end
        sampling_freq=1/(sampling_period_usec/1e6);
        voltage_trace=double(AD_count_int16)*AD_count_to_uV_factor;
        [b,a]=butter(6,filter_cutoff_frequencies/(sampling_freq/2),'bandpass'); % the Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency)
        for continuous_recording_period_i=1:length(recording_start_end_sample_indices)-1
            period_start=recording_start_end_sample_indices(continuous_recording_period_i);
            period_end=recording_start_end_sample_indices(continuous_recording_period_i+1)-1;
            filtered_voltage_traces(current_bundle_channel_i,period_start:period_end)=filtfilt(b,a,voltage_trace(period_start:period_end));
        end
    end
    clear AD_count_int16 voltage_trace
    if hard_or_adaptive_spike_threshold==1
        spike_threshold=hard_spike_threshold*ones(num_channels_on_current_bundle,1);
    elseif hard_or_adaptive_spike_threshold==2 % from Quian Quiroga et al., 2004, Neural Computation
        estimate_of_voltage_noise_std=median(abs(filtered_voltage_traces),2)/icdf('Normal',0.75,0,1);
        % Assume the fluctuations of the filtered voltage during non-spiking
        % periods (ie. the noise) is normally distributed around zero. The
        % median of the absolute values of the noise (approximated here by
        % the median of the absolute values of the entire voltage trace,
        % because spikes constitute a small proportion of the voltage
        % trace) is the 75th percentile of the noise distribution. Then,
        % this median divided by the 75th percentile of the standard normal
        % distribution (whose standard deviation is 1) is the standard
        % deviation of the noise distribution.
        spike_threshold=adaptive_spike_threshold_factor*estimate_of_voltage_noise_std;
    end
    
    if plot_voltage_traces_with_threshold
        figure(fig_threshold_check)
        num_samples_to_plot=500;
        for window_i=1:50
            clf
            samples_to_plot=(window_i-1)*num_samples_to_plot+1:window_i*num_samples_to_plot;
            times_to_plot=(1:length(samples_to_plot))*(sampling_period_usec/1000);
            for plot_channel_i=1:num_channels_on_current_bundle
                subplot(num_channels_on_current_bundle,1,plot_channel_i)
                hold on
                plot(times_to_plot,filtered_voltage_traces(plot_channel_i,samples_to_plot),'k')
                plot(times_to_plot,repmat(spike_threshold(plot_channel_i),1,length(times_to_plot)),'r--')
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
    voltage_peaks=zeros(1,num_samples_per_channel);
    for current_bundle_channel_i=1:num_channels_on_current_bundle
        [peak_heights,sample_indices_of_peaks]=findpeaks(filtered_voltage_traces(current_bundle_channel_i,:),'MinPeakHeight',spike_threshold(current_bundle_channel_i));
        peak_heights=peak_heights/spike_threshold(current_bundle_channel_i);
        higher_peaks=voltage_peaks(sample_indices_of_peaks)<peak_heights;
        voltage_peaks(sample_indices_of_peaks(higher_peaks))=peak_heights(higher_peaks);
    end
    [~,sample_indices_of_peaks]=findpeaks(voltage_peaks,'MinPeakDistance',min_separation_between_spike_peaks);
    sample_indices_of_peaks(sample_indices_of_peaks+spike_extraction_window_length(1)<1 | sample_indices_of_peaks+spike_extraction_window_length(2)>num_samples_per_channel)=[];
    clear voltage_peaks
    
    disp(['Extracting spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
    num_spikes=length(sample_indices_of_peaks);
    spike_waveforms=zeros(32,4,num_spikes);
    timestamps_usec=round(get_timestamps_for_Nlg_voltage_samples(sample_indices_of_peaks,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)); % the time stamps of all the spike peaks, rounded to integer microseconds
    for spike_i=1:num_spikes
        channel_indices_for_current_bundle=active_channels_on_current_bundle-(electrode_bundle_i-1)*channels_per_electrode_bundle;
        spike_waveforms(:,channel_indices_for_current_bundle,spike_i)=round(filtered_voltage_traces(:,sample_indices_of_peaks(spike_i)+spike_extraction_window_length(1):sample_indices_of_peaks(spike_i)+spike_extraction_window_length(2))).'; % save the waveforms of the current spike from all active channels
    end

    disp(['Saving spikes from ' lower(bundle_name) ' ' num2str(electrode_bundle_i) '...'])
    % The following works with the current version of Mat2NlxSpike (Version
    % 6.0.0) as of this writing (7/13/2016).
    NTT_file_name=[output_spike_file_name_prefix 'TT' num2str(electrode_bundle_i) '.ntt'];
    file_name_to_save=fullfile(output_folder,NTT_file_name);
    AppendToFileFlag=0;
    ExportMode=1;
    ExportModeVector=[];
    FieldSelectionFlags=[1 0 1 0 1 1]; % exporting time stamps, the clusters that spikes are assigned to, and spike waveforms
    cluster_numbers=zeros(1,length(timestamps_usec)); % zero means that the corresponding spike is not assigned to any cluster; we have to do this because Mat2NlxSpike 6.0.0 has a bug: if this input is omitted, all spikes are incorrectly assigned to clusters
    NTT_header={'######## Neuralynx Data File Header';'-ADMaxValue 500';'-ADBitVolts 1.0 1.0 1.0 1.0'};
    Mat2NlxSpike(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,timestamps_usec,cluster_numbers,spike_waveforms,NTT_header)
    % Note that if the time stamp and waveform inputs to Mat2NlxSpike
    % (6.0.0) are in certain integer formats, Mat2NlxSpike cannot correctly
    % save them, even though the .NTT file saves the waveforms and time
    % stamps in integer data formats. Thus, "timestamps_usec" and
    % "spike_waveforms" here are in the "double" format.
end
disp('Finished saving spikes from all voltage traces!')
