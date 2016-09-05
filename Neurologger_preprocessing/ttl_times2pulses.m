function [pulse_idx, pulse_time, err_pulses] = ttl_times2pulses(times,pulse_dt,correct_err,correct_end_off,correct_loop)
%%%
% Decodes spacing between TTL pulses into unique numbers
% INPUTS
% times: vector of times when TTL status changed (both on and off) in
% units of ms
%
% pulse_dt: time (ms) between pulses. Use a time less than the minimum time
% between pulses (e.g. with 5s spacing and a maximum pulse train length of
% 75ms, use approx. 4.9s; however, 4s would be sufficient)
%
% correct_err: 1 to find and correct out-of-order pulses due to e.g. 
% misreading of TTL time by NLG or rounding error. Will not correct looping
% of unique TTL files. 0 to not correct for these errors.
%
% correct_end_off: 1 to find and correct erroneous 'off' TTL at the end of
% the zero playback files. 0 to not correct for these errors.
%
% correct_loop: 1 to find and correct for looping through playback files
% (e.g. if only 1 hr of playback files are prepared, and recording lasts
% 1.5 hrs). NOTE: will only correct for ONE loop, either 0 or >1 loops will 
% throw an error. 0 to not correct for this error.
%
% NOTE: Using correct_end_off and correct_loop can be avoided by better
% constructing zero playback files (i.e. removing 'off' TTL at the end of
% files, and preparing enough total files for recording time).
%
% OUTPUTS
% pulse_idx: vector of length of the number of pulses counted. Should 
% contain the vector 1:n_pulses if decoding worked properly.
%
% pulse_time: time in ms from input 'times' when each pulse occurred.
%
% err_pulses: index of pulses in 'pulse_idx' where a pulse was out of
% order.
%
% ENCODING SCHEME: This script expects TTL pulses to be encoded in the
% following manner: every dt>pulse_dt a pulse train arrives encoding for a
% number in the sequence 1,2,3,4... always increasing by 1, but not
% necessarily starting with 1. Each pulse train is composed of a variable
% number of TTL pulses spaced by 15 ms (due to limitations imposed by NLG
% hardware). Each pulse within a pulse train is between 5 and 14 ms
% (minimum pulse length is also dictated by NLG hardware limitations). In
% order to decode these pulse trains, differences between each change in
% TTL status are calculated. Every other one of those differences are 
% extracted (i.e. we skip the 15ms spacing between pulses) and 5 is 
% subtracted from that time difference. The resulting numbers are strung
% together to form the digits of the pulse number we are decoding. 
% E.G. :
% 5s after the last pulse we see a set of TTL status change times spaced by
% 6, 15, 7, 15, and 13 ms. We look at every second number and subtract 5,
% arriving at 1,2,8 which stands for the 128th pulse in sequence. The time
% when the 6ms long pulse arrived is stored as the 'pulse_time'
% correspoding to the element of 'pulse_idx' stored as 128. 
%
% Maimon Rose 9/2/16
%%%
unique_ttls_dir = 'C:\Users\phyllo\Documents\Maimon\misc\nlg_alignment\unique_ttls\';
ttl_diffs = diff(times);
chunk_times = [times(1) times(ttl_diffs>pulse_dt) times(1)];
diffs_chunk = ttl_diffs(ttl_diffs>pulse_dt);
chunks = cell(1,length(diffs_chunk));

for chunk = 1:length(diffs_chunk)
    chunk_on = chunk_times(chunk);
    chunk_off = chunk_times(chunk+1) + 1;
    if chunk == 1
        chunks{chunk} = times((times>=chunk_on) & (times<chunk_off));
    else
        chunks{chunk} = times((times>chunk_on) & (times<chunk_off));
    end
end
pulse_time = cellfun(@(x) x(1),chunks);
display(sum(~ismember(pulse_time,times)))
chunk_diffs = cellfun(@(x) round(diff(x)),chunks,'UniformOutput',0);
pulse_idx = cellfun(@(x) str2double(regexprep(num2str((x(1:2:length(x)) - 5)),'[^\w'']','')),chunk_diffs);
pulse_idx_orig = pulse_idx;

err_pulses = intersect(find(pulse_idx-[pulse_idx(1)-1 pulse_idx(1:end-1)]~=1),... find pulses which 'stick out' from both neighboring pulses
                       find(pulse_idx-[pulse_idx(2:end) pulse_idx(end)+1]~=-1));

if correct_end_off % end of 'zero signal' TTL file may have erroneous TTL pulse at end, correct if so
    load([unique_ttls_dir 'unique_ttl_params.mat'],'n_pulse_per_chunk');
    end_offs = err_pulses(~rem(pulse_idx(err_pulses-1),n_pulse_per_chunk)); % check if erroneous TTL pulse comes at end of chunk (i.e. end of unique TTL file) by find err pulses with rem(pulse_idx,n_pulse_per_chunk) = 0
    end_offs = union(end_offs,find(isnan(pulse_idx))); 
    pulse_idx(end_offs) = []; % remove those pulses
    pulse_time(end_offs) = [];
end    
                   
err_pulses = intersect(find(pulse_idx-[pulse_idx(1)-1 pulse_idx(1:end-1)]~=1),... recalculate err_pulses after removing 'end_off' pulses
                       find(pulse_idx-[pulse_idx(2:end) pulse_idx(end)+1]~=-1));


if pulse_idx(end) - pulse_idx(end-1) ~=1
   err_pulses = [err_pulses length(pulse_idx)];
end

if correct_err
    if sum(diff(err_pulses)<2)~=0
        display([num2str(sum(diff(err_pulses)<2)) ' adjacent error pulses!']);
        try
            pulse_idx(err_pulses) = pulse_idx(err_pulses-1)+1;
        catch
            pulse_idx(err_pulses) = pulse_idx(err_pulses+1)-1;
        end
        bad_err = err_pulses(diff(err_pulses)<2);
        bad_err = [bad_err bad_err+1 bad_err-1];
        bad_err = bad_err(bad_err>1 & bad_err<=length(pulse_idx));
        pulse_idx(bad_err) = [];
    else
        try
            pulse_idx(err_pulses) = pulse_idx(err_pulses-1)+1;
        catch
            pulse_idx(err_pulses) = pulse_idx(err_pulses+1)-1;
        end
    end
end

if correct_loop
    loop_rep_idx = find(diff(pulse_idx)~=1);
    if pulse_idx(loop_rep_idx) == max(pulse_idx) && pulse_idx(loop_rep_idx+1) == min(pulse_idx)
        pulse_idx(loop_rep_idx+1:end) = pulse_idx(loop_rep_idx+1:end) + pulse_idx(loop_rep_idx);
    else
        display('no loop found!');
    end
end
figure;
hold on
plot(pulse_idx_orig,'rx');
plot(pulse_idx);
end

