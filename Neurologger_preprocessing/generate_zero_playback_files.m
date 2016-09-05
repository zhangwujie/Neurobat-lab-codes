function generate_zero_playback_files
%%%
% A script to generate so-called 'zero playback' files for synchronizing
% the avisoft audio recording system with the NLG system. This script
% generates and saves a number of large (~700 MB w/ fs = 1MHz and file
% length = 6 minutes) .WAV files which are meant to be played throug the
% avisoft player. These files have the desired TTL status encoded on their
% least significant bit (LSB) which is how the avisoft player reads out
% what TTL status it should send out. 
% ENCODING SCHEME: This script produces TTL pulses to be encoded in the
% following manner: every dt>delay a pulse train is sent out encoding for a
% number in the sequence 1,2,3,4...
% Each pulse train is composed of a variable
% number of TTL pulses spaced by 15 ms (due to limitations imposed by NLG
% hardware). Each pulse within a pulse train is between 5 and 14 ms
% (minimum pulse length is also dictated by NLG hardware limitations). In
% order to encode these pulse trains each number being encoded is broken up
% into its composite digits (e.g. 128 = [1,2,8]). We then add 5 to each of
% these values. The resulting value is the length of the TTL pulses within 
% this pulse train.  
% E.G. :
% In order to encode pulse #128 we break 128 into [1,2,8] and add 5 to
% arrive at [6,7,13]. We set the first 6 ms to TTL 'off', the next 15 ms to
% TTL 'on', the next 7 ms to TTL 'off', the next 15ms to TTL 'on', and the
% next 13ms to TTL 'off.'
% NOTE: NLG reads TTL's 'upside-down', meaning 'off' here equals 'on' for the NLG
%%%
fs = 1e6; % nominal sampling rate of avisoft player
nom_fs_offset = 93; % due to actual avisoft playback rate of 1000093
actual_fs = fs + nom_fs_offset;
[p,q] = rat(actual_fs / fs); % determine integers at which we can approximately resample data from 'fs' to 'actual_fs' such that the resulting playback files are played back at the actual fs, the playback is correct
delay = 5*fs; % interval between pulse trains
ipi = 15*fs*1e-3; % interval between pulses within individual pulse trains (set to 15ms due to NLG hardware limitations)
total_time = 1*3600*fs; % total time covered by playback files, start to finish
chunk_size = 0.1*3600*fs; % amount of time covered by one individual playback files
n_chunk = total_time/chunk_size; % number of files to be produced
n_pulse_per_chunk = chunk_size/delay - 1; % number of pulse trains in a file
n_ttl_digits = 5; % maximum number of digits in a pulse (i.e. up to 10,000 pulses)
base_ttl_length = fs*1e-3; % maximum resolution of NLG is 1ms
pulse_k = 1; % counter
for chunk = 1:n_chunk % loop through each file or 'chunk'
    wav_data = zeros(1,chunk_size,'double');
    pulse_off = []; % times when TTL status should be set of 'off' 
    for pulse = 1:n_pulse_per_chunk  % code each pulse train separately
        total_offset = delay * pulse-1; % total amount of time in this file before the pulse train about to be encoded
        ipi_offset = 0; % total amount of time taken up by the ipi's (inter pulse interval) and other pulses already within this pulse train
        d = 1; % which digit are we on
        while d <= n_ttl_digits % max pulse = 10,0000
            if pulse_k/(10^(d-1))>=1 % determine the total number of digits we need for this number
                pulse_digit_str = num2str(pulse_k); % convert pulse number to string to separate out digits
                pulse_digit = str2num(pulse_digit_str(d)); % convert the digit we are encoding now back to a number
                digit_pulse = (1:(base_ttl_length*(pulse_digit+1)+(4*base_ttl_length))) + ipi + ipi_offset; % pulse is encoded as TTL 'off' for [value of digit]*1ms + 4ms  and placed at ipi + [sum of previous ipis and pulses ] within this pulse
                pulse_off = [pulse_off,digit_pulse + total_offset]; % store this pulse along with other pulses within this file
                d = d + 1; % next digit
                ipi_offset = digit_pulse(end); % store placement for next pulse within this pulse train
            else
                d = inf; % go on to next pulse
            end
        end
        pulse_k = pulse_k + 1;
    end
    pulse_on = setdiff(1:length(wav_data),pulse_off); % set TTL to 'on' wherever it hasn't been set to 'off'
    wav_data(pulse_on) = bitset(wav_data(pulse_on),1,1); % encode the TTL status in the LSB of the data to be saved
    wav_data = int16(resample(wav_data,p,q)); % resample the data to the actual fs and convert to 16 bit integers 
    audiowrite(['unique_ttl' num2str(chunk) '.wav'],wav_data,fs); % save as a .WAV file
end
end