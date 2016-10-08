function [slope_and_intercept, mean_std_x, total_samples_by_file, first_nlg_pulse_time, first_audio_pulse_time] = align_avi_to_nlg(base_dir,ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop,wav_file_nums,session_strings)
%%
% Function to correct for clock drift between avisoft audio recordings and
% NLG neural recordings.
%
% INPUT:
% base_dir: base directory of experiment. This script expects this
% directory to contain the subfolders 'audio\ch1\' and 'nlxformat\'.
%
% ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop: see help for
% ttl_times2pulses.m
%
% wav_file_nums: vector of integers correspoding to .WAV file numbers to
% analyze.
%
% session_strings: cell of strings used to demarcate start and stop of
% time period to analyze in this script from EVENTLOG file.
%
% OUTPUT:
% slope_and_intercept: the slope and intercept of the best fit line which
% will return the estimated difference between the NLG clock and the audio
% clock (i.e. NLG time = AVI time + clock difference).
%
% mean_std_x: Centering and scaling parameters for independent variable in
% best fit line to be used in conjuction with slope_and_intercept to call
% polyval().
%
% total_samples_by_file: number of audio samples in each audio file. Used
% in order to determine time within a given audio file that is part of a
% longer recording.
%
% first_nlg_pulse_time: time (ms, in NLG time) when the first TTL pulse train
% that is used for synchronization arrived. Used to align audio and NLG
% times before scaling by estimated clock differences.
%
% first_audio_pulse_time: time (ms, in AVI time) when the first TTL pulse train
% that is used for synchronization arrived. Used to align audio and NLG
% times before scaling by estimated clock differences.
%%%
fs_wav = 250e3 + 21; % add in 21 to correct for difference between nominal avisoft clock time and actual clock time (value determined empirically)
avi_wav_bits = 16; % number of bits in each sample of avisoft data
wav2bit_factor = 2^(avi_wav_bits-1); % factor to convert .WAV data to bits readable by 'bitand' below

audio_dir = [base_dir 'audio\ch1\']; % where the audio .WAV files are stored
wav_files = dir([audio_dir '*.wav']); % all .WAV files in directory
wav_file_nums = find(cellfun(@(x) ismember(str2num(x(end-7:end-4)),wav_file_nums),{wav_files.name})); % extract only requested .WAV files
audio_time_din = [];
total_samples = 0;
total_samples_by_file = zeros(1,length(wav_files));

save_options_parameters_CD_figure = 1;

for w = 1:max(wav_file_nums) % run through all requested .WAV files and extract audio data and TTL status at each sample
    if ismember(w,wav_file_nums)
        data = audioread([audio_dir wav_files(w).name]); % load audio data
        ttl_status = bitand(data*wav2bit_factor + wav2bit_factor,1); % read TTL status off least significant bit of data
        audio_time_din = [audio_time_din (1e3*(total_samples + find(sign(diff(ttl_status))~=0)')/fs_wav)];
        total_samples_by_file(w) = length(data);
        total_samples = total_samples + total_samples_by_file(w);
    else
        audio_info_struct = audioinfo([audio_dir wav_files(w).name]);
        total_samples_by_file(w) = audio_info_struct.TotalSamples;
        total_samples = total_samples + total_samples_by_file(w);
    end
end

[audio_pulses, audio_pulse_times] = ttl_times2pulses(audio_time_din,ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop); % extract TTL pulses and time

%%

eventfile = [base_dir 'nlxformat\EVENTS.mat']; % load file with TTL status info
csc_file = [base_dir 'nlxformat\CSC0.mat'];
load(eventfile);
load(csc_file,'indices_of_first_samples','timestamps_of_first_samples_usec','sampling_period_usec','indices_of_first_samples')

session_start = event_timestamps_usec(find(cellfun(@(x) ~isempty(strfind(x,session_strings{1})),event_types_and_details),1,'last')); % find when to start analyzing
session_end = event_timestamps_usec(find(cellfun(@(x) ~isempty(strfind(x,session_strings{2})),event_types_and_details),1,'last')); % find when to stop analyzing
if isempty(session_start)
    display('couldn''t find start session string in event file, choose index of events to use as session start');
    keyboard;
    session_start_idx = input('input index into variable event_types_and_details');
    session_start = event_timestamps_usec(session_start_idx);
end
if isempty(session_end)
    display('couldn''t find end session string in event file, choose index of events to use as session start');
    keyboard;
    session_end_idx = input('input index into variable event_types_and_details');
    session_end = event_timestamps_usec(session_end_idx);
end

% n_nlg_files = length(timestamps_of_first_samples_usec); % total number of NLG files created
% nlg_indices = ((session_start - timestamps_of_first_samples_usec(1) - n_nlg_files*1e3)/sampling_period_usec):((session_end - timestamps_of_first_samples_usec(1) + n_nlg_files*1e3)/sampling_period_usec); %
% nlg_timestamps_usec = get_timestamps_for_Nlg_voltage_samples(nlg_indices,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec); % build vector of timestamps (usec) for all NLg samples

% extract only relevant TTL status changes
event_types_and_details = event_types_and_details((event_timestamps_usec >= session_start(end)) & (event_timestamps_usec <= session_end(end)));
event_timestamps_usec = event_timestamps_usec((event_timestamps_usec >= session_start(end)) & (event_timestamps_usec <= session_end(end)));

din = cellfun(@(x) ~isempty(strfind(x,'Digital in')),event_types_and_details); % extract which lines in EVENTS correspond to TTL status changes
nlg_time_din = 1e-3*event_timestamps_usec(din)'; % find times (ms) when TTL status changes
[nlg_pulse, nlg_pulse_times] = ttl_times2pulses(nlg_time_din,ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop); % extract TTL pulses and time

%% synchronize audio --> NLG
[~, shared_pulse_nlg_idx, shared_pulse_audio_idx] = intersect(nlg_pulse,audio_pulses); % determine which pulses are on both the NLG and avisoft recordings

% extract only shared pulses
shared_nlg_pulse_times = nlg_pulse_times(shared_pulse_nlg_idx);
shared_audio_pulse_times = audio_pulse_times(shared_pulse_audio_idx);

first_nlg_pulse_time = shared_nlg_pulse_times(1);
first_audio_pulse_time = shared_audio_pulse_times(1);

% nlg_t = 1e-3*nlg_timestamps_usec; % convert timestamps to ms
% nlg_t_idx = (nlg_t>=shared_nlg_pulse_times(1)) & (nlg_t<=shared_nlg_pulse_times(end)); % extract timestamps only from time when NLG and avisoft are sharing pulses
% nlg_t = nlg_t(nlg_t_idx);
%
% audio_t = (1e3*((0:length(audio_ttl))/fs_wav)); % convert time to msec
% audio_t_idx = (audio_t>=shared_audio_pulse_times(1)) & (audio_t<=shared_audio_pulse_times(end)); % extract timestamps only from time when NLG and avisoft are sharing pulses
% audio_t = audio_t(audio_t_idx);

clock_differences_at_pulses = (shared_nlg_pulse_times - first_nlg_pulse_time) - (shared_audio_pulse_times - first_audio_pulse_time); % determine difference between NLG and avisoft timestamps when pulses arrived

[slope_and_intercept,~,mean_std_x]=polyfit(shared_audio_pulse_times-first_audio_pulse_time,clock_differences_at_pulses,1); % fit line to all clock differences with respect to avisoft timestamps

estimated_clock_diff =  polyval(slope_and_intercept,shared_audio_pulse_times - first_audio_pulse_time,[],mean_std_x);
figure
subplot(2,1,1)
hold on
plot(clock_differences_at_pulses,clock_differences_at_pulses-estimated_clock_diff,'x')
xlabel('difference between NLG clock and avisoft clock');
ylabel('deviation of estimated clock differences from real clock differences');
subplot(2,1,2)
hold on
plot(shared_audio_pulse_times-first_audio_pulse_time,clock_differences_at_pulses,'.-');
plot(shared_audio_pulse_times-first_audio_pulse_time,estimated_clock_diff,'r');
xlabel('Incoming Audio Pulse Times')
ylabel('Difference between NLG clock and avisoft clock');
legend('real clock difference','estimated clock difference');

if save_options_parameters_CD_figure
    saveas(gcf,fullfile(audio_dir,'CD_correction_avisoft_nlg.fig'))
end

% corrected_audio_t = audio_t + estimated_clock_differences; % add those differences to avisoft timestamps
%
% % reference each set of timestamps to the first TTL pulse received on each
% % system
% corrected_audio_t  = corrected_audio_t - corrected_audio_t(1);
% nlg_t = nlg_t - nlg_t(1);

end



