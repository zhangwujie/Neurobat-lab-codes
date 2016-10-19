function corr_t = avi2nlg_time(audio2nlg,t)
% Helper function to convert from Avisoft time to NLG time.
% INPUT:
%
% audio2nlg: Structure of outputs from align_avi_to_nlg with fields:
%   'shared_nlg_pulse_times','shared_audio_pulse_times','total_samples_by_file','first_audio_pulse_time','first_nlg_pulse_time'
%
% t: Avisoft time in ms
%
% OUTPUT:
%
% corr_t: NLG time in ms

t = t - audio2nlg.first_audio_pulse_time;
clock_differences_at_pulses = (audio2nlg.shared_nlg_pulse_times - audio2nlg.first_nlg_pulse_time) - (audio2nlg.shared_audio_pulse_times - audio2nlg.first_audio_pulse_time);
estimated_clock_differences = interp1(audio2nlg.shared_audio_pulse_times,clock_differences_at_pulses,t);
corr_t = t + estimated_clock_differences;
end
