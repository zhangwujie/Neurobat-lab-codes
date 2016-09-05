function [nlg_t, corrected_audio_t, audio_data, audio_t_idx, nlg_csc, nlg_t_idx] = align_avi_to_nlg(base_dir,ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop,wav_file_nums,session_strings)
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
% nlg_t: vector of timestamps for each sample collected by NLG in ms.
%
% corrected_audio_t: vector of timestamps for each sample collected by
% avisoft in ms, corrected for drift between NGL and avisoft.
%
% audio_data: raw audio data from all audio files considered concatenated
% into one vector within bounds imposed by 'session_strings' and first TTL
% pulse.
%
% audio_t_idx: index into concatenation of all .WAV files considered which
% produces 'audio_data'
%
% nlg_csc: NLG data from full recording session within bounds imposed by 
% 'session_strings' and first TTL pulse.
%
% nlg_t_idx: index into NLG data which  produces 'nlg_data'
%%%
fs_wav = 250e3 + 21; % add in 21 to correct for difference between nominal avisoft clock time and actual clock time (value determined empirically)
avi_wav_bits = 16; % number of bits in each sample of avisoft data
wav2bit_factor = 2^(avi_wav_bits-1); % factor to convert .WAV data to bits readable by 'bitand' below

audio_dir = [base_dir 'audio\ch1\']; % where the audio .WAV files are stored
wav_files = dir([audio_dir '*.wav']); % all .WAV files in directory
wav_file_nums = find(cellfun(@(x) ismember(str2num(x(end-7:end-4)),wav_file_nums),{wav_files.name})); % extract only requested .WAV files
audio_data = [];
audio_ttl = [];
total_samples = 0;

for w = 1:max(wav_file_nums) % run through all requested .WAV files and extract audio data and TTL status at each sample
    if ismember(w,wav_file_nums)
        data = audioread([audio_dir wav_files(w).name]); % load audio data
        ttl_status = bitand(data*wav2bit_factor + wav2bit_factor,1); % read TTL status off least significant bit of data
        audio_ttl = [audio_ttl;ttl_status]; % store TTL status
        audio_data = [audio_data;data]; % store audio data
    else
        audio_info_struct = audioinfo([audio_dir wav_files(w).name]);
        total_samples = total_samples + audio_info_struct.TotalSamples;
        ttl_status = ones(total_samples,1);
        audio_ttl = [audio_ttl;ttl_status];
        audio_data = [audio_data; zeros(total_samples,1)];
    end
end

audio_time_din = 1e3*(find(sign(diff(audio_ttl))~=0)')/fs_wav; % find samples where TTL status changes and convert from samples to ms
[audio_pulses, audio_pulse_times] = ttl_times2pulses(audio_time_din,ttl_pulse_dt,corr_pulse_err,correct_end_off,correct_loop); % extract TTL pulses and time

%%

eventfile = [base_dir 'nlxformat\EVENTS.mat']; % load file with TTL status info
CSC_file = [base_dir 'nlxformat\CSC0.mat']; % load file with NLG data in order to create timestamps
load(eventfile);
load(CSC_file);
nlg_timestamps_usec = get_timestamps_for_Nlg_voltage_samples(1:length(AD_count_int16),indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec); % build vector of timestamps (usec) for all NLg samples
session_start = event_timestamps_usec(cellfun(@(x) strcmp(x,['Free text. ' session_strings{1}]),event_types_and_details)); % find when to start analyzing
session_end = event_timestamps_usec(cellfun(@(x) strcmp(x,['Free text. ' session_strings{2}]),event_types_and_details)); % find when to stop analyzing

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

nlg_t = 1e-3*nlg_timestamps_usec; % convert timestamps to ms
nlg_t_idx = (nlg_t>=shared_nlg_pulse_times(1)) & (nlg_t<=shared_nlg_pulse_times(end)); % extract timestamps only from time when NLG and avisoft are sharing pulses 
nlg_t = nlg_t(nlg_t_idx);

audio_t = (1e3*((0:length(audio_ttl))/fs_wav)); % convert time to msec
audio_t_idx = (audio_t>=shared_audio_pulse_times(1)) & (audio_t<=shared_audio_pulse_times(end)); % extract timestamps only from time when NLG and avisoft are sharing pulses 
audio_t = audio_t(audio_t_idx);

clock_differences_at_pulses = (shared_nlg_pulse_times - shared_nlg_pulse_times(1)) - (shared_audio_pulse_times - shared_audio_pulse_times(1)); % determine difference between NLG and avisoft timestamps when pulses arrived

[slope_and_intercept,~,mean_std_x]=polyfit(shared_audio_pulse_times,clock_differences_at_pulses,1); % fit line to all clock differences with respect to avisoft timestamps
estimated_clock_differences = polyval(slope_and_intercept,audio_t,[],mean_std_x); % evaluate that line at all avisoft timestamps to estimate difference between NLG and avisoft clocks at every sample
corrected_audio_t = audio_t + estimated_clock_differences; % add those differences to avisoft timestamps

% reference each set of timestamps to the first TTL pulse received on each
% system
corrected_audio_t  = corrected_audio_t - corrected_audio_t(1); 
nlg_t = nlg_t - nlg_t(1);

nlg_csc = double(AD_count_int16(nlg_t_idx)*AD_count_to_uV_factor);
audio_data = audio_data(audio_t_idx);

end



