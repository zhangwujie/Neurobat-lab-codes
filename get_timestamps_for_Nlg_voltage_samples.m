% When using extract_Nlg_data.m to extract the voltage trace data from
% Neurologger .DAT file and saving as .mat files, we do not save a time
% stamp for each voltage sample. Instead, we save all information needed to
% calculate the time stamps in the same .mat file as the voltage data, and
% use this function to calculate the time stamps for the samples we're
% interested in.
% 7/12/2016, Wujie Zhang

function timestamps_usec=get_timestamps_for_Nlg_voltage_samples(sample_indices,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)
timestamps_usec=nan(size(sample_indices));
for i=1:numel(sample_indices)
    file_index=find(sample_indices(i)>=indices_of_first_samples,1,'last');
    periods_from_first_sample_in_file=sample_indices(i)-indices_of_first_samples(file_index);
    timestamps_usec(i)=timestamps_of_first_samples_usec(file_index)+periods_from_first_sample_in_file*sampling_period_usec;
end