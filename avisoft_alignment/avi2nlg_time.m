function corr_t = avi2nlg_time(audio2nlg,t,method)
% Helper function to convert from Avisoft time to NLG time.
% INPUT:
%
% audio2nlg: Structure of outputs from align_avi_to_nlg with fields:
%   'shared_nlg_pulse_times','shared_audio_pulse_times','total_samples_by_file','first_audio_pulse_time','first_nlg_pulse_time'
%
% t: Avisoft time in ms, counting from the first sample of the first audio file
%
% method: 1 means fitting a single line over all points, 2 means
% interpolating between consecutive points
%
% OUTPUT:
%
% corr_t: NLG time in ms, counting from the time of the first TTL chunk,
% which is audio2nlg.first_nlg_pulse_time
% 
% Maimon Rose
% Last updated: 11/1/2016, Wujie Zhang

t = t - audio2nlg.first_audio_pulse_time;
clock_differences_at_pulses = (audio2nlg.shared_nlg_pulse_times - audio2nlg.first_nlg_pulse_time) - (audio2nlg.shared_audio_pulse_times - audio2nlg.first_audio_pulse_time);

if method==1
    [slope_and_intercept,~,mean_std_x]=polyfit(audio2nlg.shared_audio_pulse_times,clock_differences_at_pulses,1);
    estimated_clock_differences=polyval(slope_and_intercept,t,[],mean_std_x);
elseif method==2
    estimated_clock_differences = interp1(audio2nlg.shared_audio_pulse_times,clock_differences_at_pulses,t,'linear','extrap');
end

corr_t = t + estimated_clock_differences;
end
