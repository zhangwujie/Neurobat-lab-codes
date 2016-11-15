% Given a folder containing Neurologger (Nlg) data from one recording
% session, extract the event and voltage data, and save in either the
% MATLAB .mat format, or in the Neuralynx (Nlx) .nev format (for the event
% file) and .ncs format (for the raw voltage traces represented by
% analog-to-digital counts [AD counts] in bits) in another folder.
% -The MATLAB format is recommended over the Nlx formats, because it's
% simpler, easier to work with, has less bugs, more thoroughly tested, and
% has much smaller file sizes.
% -The Neurologger event file contains time stamps from the Neurologger and
% its transceiver, with clock drift between the two -- we correct the clock
% drift so that all time stamps are those of the Neurologger clock.
% -Can choose whether or not to save event and voltage data files
% -Before running, need to convert the Neurologger .NLE event file to .xlsx
% format
% 7/2/2016, Wujie Zhang
last_code_update='9/6/2016, Wujie Zhang'; % identifies the version of the code
%%
% Input and ouput paths, options, and paramters
Nlg_folders={'F:\Wujie\Data\bat selection recording 1\20161018 cf\nlgformat\'}; % each cell is a folder where Nlg voltage data files and event file are stored
output_folders_for_each_Nlg_folder={'F:\Wujie\Data\bat selection recording 1\20161018 cf\neural\'}; % each cell is a folder where the outputs of the code (voltage and event data in MATLAB or Nlx formats) will be saved, corresponding to one of the folders in Nlg_folders

inactive_channels_for_each_Nlg_folder={[2:16]}; % each cell is a vector containing numbers between 1 and the number of channels, indicating the disabled channels; enter an empty vector if all channels are active; each cell corresponds to one of the folders in Nlg_folders
reference_channel_for_each_Nlg_folder={[]}; % each cell is a number between 1 and the number of channels, indicating a reference channel whose AD counts (equivalently, voltages) will be subtracted from those of all other channels during the processing here; enter an empty vector to not subtract any reference channel here; each cell corresponds to one of the folders in Nlg_folders

save_in_mat_or_Nlx_format=1; % 1: save in .mat format; 2: save in Nlx .nev and .ncs formats
save_event_file=1; % 0: don't save any event file; 1: save one event file containing all events; 2: save one event file containing all events and one event file for each event type
save_voltage_AD_count_files=1; % 0: don't save voltage data files; 1: save
save_options_parameters_CD_figure=1; % 1: save the options and paramters for the current run of the code in a .mat file, and the figure showing the clock difference correction; 0: don't save
global_or_local_clock_difference_estimation=1; % to estimate unreported clock differences, 1: fit a line using all reported clock differences in a recording interval where clocks are not synchronized; 2: interpolate between consecutive time points with reported clock differences
timestamps_from_eventlog_or_sampling_period=2; % for the time stamps of the first sample of each Nlg .DAT file, 1: take them from the event log; 2: calculate them by counting the number of samples since recording start, and multiplying by the sampling period; the second method might be more accurate, because the event log time stamps are rounded to integer milliseconds, while the sampling period is tens of microseconds
%%
% Parameters specific to the Neurologger used
Nlg_file_name_letters='NEUR'; % the first four letters of the Nlg .DAT file names; eg. 'NEUR' for file names like "NEUR_003.DAT"
num_channels=16; % total number of recording channels, including inactive ones
AD_count_for_zero_voltage=2048; % the AD count that represents zero volt, which will be subtracted from the recorded AD counts during the processing here
AD_count_to_uV_factor=3.3; % voltage in microvolt = (raw AD count - AD_count_for_zero_voltage) * AD_count_to_uV_factor; this code does not convert the AD counts to voltages--this conversion factor is saved along with the AD counts in the same .mat files, to be used by detect_spikes_from_raw_voltage_trace.m; alternatively, leave as [] to read this from the Nlg event file header
sampling_period_sec=(512/15)/1e6; % specify here the exact sampling period in s for the signal in a given recording channel, or leave it as [] to read the sampling period from the Nlg event file header; this is the number of channels times the sampling period of the AD converter
% For the Neurolog-16 that we use as of this writing (7/4/2016),
% 512/15 = 34.13333... us is exactly 16 times the sampling period of the AD
% converter. The sampling period from the header (34.133328 us) is slightly
% inaccurate (would result in a difference of ~1 ms over ~2 hours compared
% to using 512/15 us).
samples_per_channel_per_file=524288; % specify here the number of samples per channel per Nlg .DAT file; alternatively, leave as [] to automatically calculate this from the number of samples in the first file loaded and the number of channels, num_channels
% For the Neurolog-16 that we use as of this writing (7/6/2016), each file
% contain 2^24 bytes, and each sample takes 2 bytes, so the samples per
% channel per file is 2^24 bytes / (2 bytes * 16 channels) = 524288
Nlg_data_type='uint16'; % the data type of the AD counts in the Nlg .DAT files
% For the Neurolog-16 that we use as of this writing (7/11/2016), the data
% type is 16-bit unsigned integer.
Nlg_unwritten_data_value=65535; % specify here the default value of Nlg voltage data samples when they haven't been written with actual recorded data
% When recording is stopped (at the end of a recording session or in the
% middle), the latest recorded data are written to the latest .DAT file,
% which partially fills it and leaves all samples in the unfilled portion
% of the file at their default value, which is 65535 (the 16-bit unsigned
% integer with all bits being ones) in the version of Neurolog-16 used as
% of this writing (7/11/2016), and is 0 in some other versions of
% Neurologgers.
old_clock_difference_bug_correction=0; % if processing data recorded when the Nlg had a bug in reporting the clock difference (ie. bat #95, "Primus", from early 2016), set this to 1 to correct the error
%%
for Nlg_folder_i=1:length(Nlg_folders) % for each of the Nlg folders
    Nlg_folder=Nlg_folders{Nlg_folder_i};
    output_folder=output_folders_for_each_Nlg_folder{Nlg_folder_i};
    inactive_channels=inactive_channels_for_each_Nlg_folder{Nlg_folder_i};
    reference_channel=reference_channel_for_each_Nlg_folder{Nlg_folder_i};
    disp(['Processing the Nlg data in "' Nlg_folder '"...'])
    if save_event_file || save_voltage_AD_count_files || save_options_parameters_CD_figure % if we're saving anything
        if ~exist(output_folder,'dir'); % make the output folder if it doesn't already exist
            mkdir(output_folder);
        end
    end
    %%
    % read Nlg event file; this code assumes that the columns of the Nlg event
    % file are, in order, "Event Number", "Time Stamp", "Time (ms from
    % midnight)", "Time Source", "Event Type", and "Details" (7/20/2016)
    Nlg_event_file_xlsx_name=fullfile(Nlg_folder,'EVENTLOG.xlsx'); % file path of the Nlg event file in .xlsx format
    [~,~,imported_Nlg_event_xlsx]=xlsread(Nlg_event_file_xlsx_name); % importing the entire excel sheet as a cell array, importing numbers as numbers, words as strings, and empty cells as NaNs
    string_or_number=zeros(size(imported_Nlg_event_xlsx,1),1);
    for row_i=1:size(imported_Nlg_event_xlsx,1) % looping over the rows of the spreadsheet (which is in general faster than using "cellfun"), to find out whether the entry in the first column is string or number
        if ischar(imported_Nlg_event_xlsx{row_i,1})
            string_or_number(row_i)=1;
        elseif isnumeric(imported_Nlg_event_xlsx{row_i,1}) && ~isnan(imported_Nlg_event_xlsx{row_i,1})
            string_or_number(row_i)=2;
        else
            disp(['Row ' num2str(row_i) ', column 1 in excel sheet is neither number nor string.'])
        end
        if isnan(imported_Nlg_event_xlsx{row_i,6}) % when "Details" for a row is empty, MATLAB imports it as Nan
            imported_Nlg_event_xlsx{row_i,6}=''; % change the NaN to an empty string
        end
    end
    first_numeric_row=find(string_or_number==2,1); % the first few rows of the spreadsheet contain the header and the labels for the columns, and the next rows contain the events, which have numbers in the first column; this line finds the first of those rows with numbers
    
    continued_event_indices=find(string_or_number==1);
    continued_event_indices=continued_event_indices(continued_event_indices>first_numeric_row); % find the rows that are "...Continued" from the previous rows, which have the string "-----" in the first column
    for continued_i=continued_event_indices % for each of the "...Continued" events, add its "Details" to that of the event from which it contnues
        continued_from_i=find(string_or_number(1:continued_i)==2,1,'last'); % the last row before the current "...Continued" row that has a number, ie. the event from which the current "...Continued" row continues
        imported_Nlg_event_xlsx{continued_from_i,6}=[imported_Nlg_event_xlsx{continued_from_i,6} imported_Nlg_event_xlsx{continued_i,6}];
    end
    imported_Nlg_event_xlsx(continued_event_indices,:)=[]; % delete all the "...Continued" rows, since their "Details" have already been copied
    
    if isempty(sampling_period_sec) % if the sampling period wasn't manually entered above, read it from the header
        Nlg_header_rows=1:first_numeric_row-2; % the row before the first numeric row is the column labels (eg. 'Time Stamp', 'Event Type', etc.), and the rows before that are the header
        Nlg_header=[imported_Nlg_event_xlsx{Nlg_header_rows,1}];
        string_before='ADC Period = ';
        string_after='us';
        string_start_index=strfind(Nlg_header,string_before)+length(string_before);
        string_end_index=strfind(Nlg_header,string_after)-1;
        sampling_period_sec=str2double(Nlg_header(string_start_index:string_end_index))/1e6;
    end
    if isempty(AD_count_to_uV_factor) % if the conversion factor from AD counts to voltages wasn't manually entered above, read it from the header
        Nlg_header_rows=1:first_numeric_row-2; % the row before the first numeric row is the column labels (eg. 'Time Stamp', 'Event Type', etc.), and the rows before that are the header
        Nlg_header=[imported_Nlg_event_xlsx{Nlg_header_rows,1}];
        string_before='ADC Resolution = ';
        string_after='uV';
        string_start_index=strfind(Nlg_header,string_before)+length(string_before);
        string_end_index=strfind(Nlg_header,string_after)-1;
        AD_count_to_uV_factor=str2double(Nlg_header(string_start_index:string_end_index));
    end
    
    event_timestamps_usec=cell2mat(imported_Nlg_event_xlsx(first_numeric_row:end,3))*1000; % time stamps, converted from ms to us
    event_timestamps_source=imported_Nlg_event_xlsx(first_numeric_row:end,4); % source of time stamps: logger or transceiver
    event_types=imported_Nlg_event_xlsx(first_numeric_row:end,5);
    event_types_and_details=strcat(imported_Nlg_event_xlsx(first_numeric_row:end,5),{'. '},imported_Nlg_event_xlsx(first_numeric_row:end,6)); % event type and details are concatenated together, with a period and space between them, eg. "File started" and "File index: 001" and concatenated into "File started. File index: 001"
    %%
    % Change all transceiver time stamps to the corresponding times of the
    % logger clock (the transceiver clock and logger clock run at different
    % speeds); the clock difference is defined as logger time minus transceiver
    % time
    string_before_CD='CD=';
    pre_CD_positions=strfind(event_types_and_details,string_before_CD);
    PC_comments_indices=find(~cellfun(@isempty,pre_CD_positions)); % find all the "PC-generated comment" events
    clock_differences_sec=nan(length(PC_comments_indices),1);
    logger_times=nan(length(PC_comments_indices),1);
    for PC_comment_i=1:length(PC_comments_indices)
        current_row_i=PC_comments_indices(PC_comment_i);
        string_after_CD=event_types_and_details{current_row_i}(pre_CD_positions{current_row_i}+length(string_before_CD):end);
        position_CD_end=find(string_after_CD==' ',1)-1; % find the empty space after the clock difference value from eg. "...CD=-0.001000 RR=0..."
        clock_differences_sec(PC_comment_i)=str2double(string_after_CD(1:position_CD_end)); % clock differences are in s
        logger_times(PC_comment_i)=event_timestamps_usec(current_row_i); % the time of the logger clock when the clock difference was reported
    end
    if old_clock_difference_bug_correction % only use this for data collected when transciever had error, ie. Mindy's data from bat #95, "Primus"
        clock_differences_sec=clock_differences_sec*16+0.03;
    end
    clock_differences_usec=clock_differences_sec*1e6; % convert from s to us
    transceiver_times=logger_times-clock_differences_usec; % the times of the transceiver clock when the clock differences were reported
    
    indices_bordering_unsynchronized_intervals=[1; find(ismember(event_types_and_details,'Clocks synchronized. ')); length(event_types_and_details)]; % the events when the two clocks are synchronzied, and the first and last event
    indices_events_with_transceiver_time=find(ismember(event_timestamps_source,'Transceiver')); % find the events that were originally logged with transceiver time stamps
    estimated_clock_differences=nan(length(indices_events_with_transceiver_time),1);
    
    for unsynchronized_interval_i=1:length(indices_bordering_unsynchronized_intervals)-1 % for each of the intervals between consecutive "Clocks synchronized" events
        logical_indices_transc_t_in_interval=indices_events_with_transceiver_time>indices_bordering_unsynchronized_intervals(unsynchronized_interval_i) & indices_events_with_transceiver_time<indices_bordering_unsynchronized_intervals(unsynchronized_interval_i+1); % all the events that were originally logged with transceiver time stamps in the current interval
        if ~any(logical_indices_transc_t_in_interval) % if there is no events in the current interval whose clock difference need to be estimated
            continue
        end
        logical_indices_PC_comments_in_interval=PC_comments_indices>indices_bordering_unsynchronized_intervals(unsynchronized_interval_i) & PC_comments_indices<indices_bordering_unsynchronized_intervals(unsynchronized_interval_i+1); % all the events in the current interval when clock difference was reported
        transceiver_times_in_interval=transceiver_times(logical_indices_PC_comments_in_interval);
        clock_differences_usec_in_interval=clock_differences_usec(logical_indices_PC_comments_in_interval);
        if length(clock_differences_usec_in_interval)==1 % if there is only one reported clock difference in the current interval, then use that for all unreported clock differences
            estimated_clock_differences(logical_indices_transc_t_in_interval)=clock_differences_usec_in_interval;
        else
            % Estimate the clock differences at all the time points in the
            % current interval that were originally logged with the
            % transceiver clock
            if global_or_local_clock_difference_estimation==1
                [slope_and_intercept,~,mean_std_x]=polyfit(transceiver_times_in_interval,clock_differences_usec_in_interval,1);
                estimated_clock_differences(logical_indices_transc_t_in_interval)=polyval(slope_and_intercept,event_timestamps_usec(indices_events_with_transceiver_time(logical_indices_transc_t_in_interval)),[],mean_std_x);
                % Estimate by fitting a line over the clock differences of
                % all the PC-generated comments in the current interval
            elseif global_or_local_clock_difference_estimation==2
                estimated_clock_differences(logical_indices_transc_t_in_interval)=interp1(transceiver_times_in_interval,clock_differences_usec_in_interval,event_timestamps_usec(indices_events_with_transceiver_time(logical_indices_transc_t_in_interval)),'linear','extrap');
                % Estimate by linearly interpolating between the clock
                % differences of each pair of consecutive PC-generated
                % comments; also extrapolates for time points before the
                % first PC-generated comment or after the last
            end
        end
    end
    if any(isnan(estimated_clock_differences))
        disp('Unable to estimate clock differences at some of the events...')
    end
    event_timestamps_usec(indices_events_with_transceiver_time)=event_timestamps_usec(indices_events_with_transceiver_time)+estimated_clock_differences; % convert the time stamps that were originally transceiver times to logger times
    
    event_timestamps_usec=round(event_timestamps_usec); % round all time stamps to integer microseconds
    
    figure % plot the result of clock difference correction to check for mistakes
    hold on
    plot((event_timestamps_usec(indices_bordering_unsynchronized_intervals(2:end-1))-event_timestamps_usec(1))/(1e6*60),zeros(length(indices_bordering_unsynchronized_intervals(2:end-1)),1),'k*') % the events when the two clocks are synchronized
    plot((logger_times-event_timestamps_usec(1))/(1e6*60),clock_differences_usec/1e3,'ro') % logger times vs. the reported clock differences that were used for estimation
    plot((event_timestamps_usec(indices_events_with_transceiver_time)-event_timestamps_usec(1))/(1e6*60),estimated_clock_differences/1e3,'b.') % logger times vs. the estimated clock differences for all the time stamps that were originally transceiver times
    legend('Clock synchronization events','Recorded clock differences','Estimated clock differences','Location','best')
    % Even right after clock synchronization, the reported clock difference
    % may not be zero, because there is a finite uncertainty in the
    % reported clock difference
    hold on
    ylabel('Logger time - transceiver time (ms)')
    xlabel('Logger time (minutes)')
    if save_options_parameters_CD_figure
        saveas(gcf,fullfile(output_folder,'CD_correction.fig'))
    end
    %%
    % save event files in either MATLAB or Nlx format
    if save_event_file
        if save_in_mat_or_Nlx_format==1 % if saving in MATLAB format
            file_name_to_save=fullfile(output_folder,'EVENTS.mat');
            clear variables_to_save
            variables_to_save.event_timestamps_usec=event_timestamps_usec;
            variables_to_save.event_types_and_details=event_types_and_details;
            save(file_name_to_save,'-struct','variables_to_save')
            
            if save_event_file==2 % if also saving sperate event files for each event type
                event_type_list=unique(event_types); % all the event types
                clear variables_to_save
                for event_type_i=1:length(event_type_list)
                    current_event_type=event_type_list{event_type_i};
                    indices_of_events_of_this_type=ismember(event_types,current_event_type); % find all the events of this type
                    file_name_to_save=fullfile(output_folder,['EVENTS_' current_event_type '.mat']);
                    variables_to_save.event_timestamps_usec=event_timestamps_usec(indices_of_events_of_this_type);
                    variables_to_save.event_types_and_details=event_types_and_details(indices_of_events_of_this_type);
                    save(file_name_to_save,'-struct','variables_to_save')
                end
            end
        elseif save_in_mat_or_Nlx_format==2 % if saving in Nlx format
            file_name_to_save=fullfile(output_folder,'EVENTS.nev');
            
            % The following works with the current version of Mat2NlxEV
            % (Version 6.0.0) as of this writing (7/5/2016).
            AppendToFileFlag=0; % save new file, or overwrite existing file, but do not append to existing file
            ExportMode=1; % export all events
            ExportModeVector=[]; % don't need this when exporting all events
            FieldSelectionFlags=[1 0 0 0 1 0]; % whether to export: Timestamps, Event IDs, TTLs, Extras, Event Strings, Header; the Neurologger TTLs are treated as regular events, so don't export TTLs here
            Mat2NlxEV(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,event_timestamps_usec.',event_types_and_details);
            % For Mat2NlxEV 6.0.0, the "Timestamps" input needs to be 1 X N,
            % and the "Event Strings" input needs to be N X 1 (unlike what its
            % documentation states); note that the Nlx .nev format only stores
            % time in integer microseconds, and Mat2NlxEV rounds non-integer
            % times towards zero (instead of towards the nearest integer),
            % so make sure all time stamps are integers before exporting
            % (as they are in this code)
            
            if save_event_file==2 % if also saving sperate event files for each event type
                event_type_list=unique(event_types); % all the event types
                for event_type_i=1:length(event_type_list)
                    current_event_type=event_type_list{event_type_i};
                    indices_of_events_of_this_type=ismember(event_types,current_event_type); % find all the events of this type
                    file_name_to_save=fullfile(output_folder,['EVENTS_' current_event_type '.nev']);
                    Mat2NlxEV(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,event_timestamps_usec(indices_of_events_of_this_type).',event_types_and_details(indices_of_events_of_this_type));
                end
            end
        end
    end
    %%
    % Save the AD counts (in bits) encoding the raw voltage traces
    if save_voltage_AD_count_files
        indices_missing_Nlg_files=[];
        ADC_sampling_period_usec=(sampling_period_sec*1e6)/num_channels; % the Nlg AD converter samples each channel sequentially, so its sampling period is the sampling period of a single channel divided by the number of channels
        active_channels=1:num_channels;
        active_channels(inactive_channels)=[]; % take out the inactive channels
        num_active_channels=length(active_channels);
        
        indices_file_start_events=find(ismember(event_types,'File started')); % find the 'File started' events
        file_start_timestamps_usec=event_timestamps_usec(indices_file_start_events);
        file_start_details=event_types_and_details(indices_file_start_events); % these are eg. "File started. File index: 000"
        num_Nlg_files=length(file_start_timestamps_usec);
        
        indices_stop_recording_events=find(ismember(event_types_and_details,'Mode change. Stopped recording'));
        find_partially_filled_Nlg_files=zeros(length(indices_stop_recording_events),1);
        for stop_recording_i=1:length(indices_stop_recording_events) % for each of the "Stopped recording" events
            current_stop_recording_index=indices_stop_recording_events(stop_recording_i);
            find_partially_filled_Nlg_files(stop_recording_i)=max(indices_file_start_events(indices_file_start_events<current_stop_recording_index)); % the last "File started" event before the "Stopped recording" event
        end
        find_partially_filled_Nlg_files=[find_partially_filled_Nlg_files; indices_file_start_events(end)]; %#ok<AGROW> % add in the last .DAT file, because if the battery runs out or the logger is turned off before the "stop recording" button is pressed, this last file is not followed by a "Stopped recording" event
        find_partially_filled_Nlg_files=unique(find_partially_filled_Nlg_files); % find unique indices, because the "Stopped recording" event sometimes happen twice in a row without any recording in-between
        found_event_strings=cell2mat(event_types_and_details(find_partially_filled_Nlg_files)); % these are eg. "File started. File index: 000"
        partially_filled_Nlg_files=find(ismember(file_start_details,found_event_strings)); % the indices of the partially filled files; note these indices start from 1, so that 1 is Nlg .DAT file 000, 2 is Nlg .DAT file 001, etc.
        
        first_file_loaded=1;
        first_file_after_recording_start=1; % as the code processes each Nlg .DAT file in sequence, this variable will be updated to reflect the index of the first file during the current period of continuous recording (eg. if recording was stopped and restarted during one experimental session, this would be the first file after recording was restarted)
        file_start_timestamp_discrepancies_ms=nan(num_Nlg_files,1); % for each .DAT file, the file start time stamp read from the event log minus that calculated by number of samples from recording start, after rounding both time stamps to integer ms
        for Nlg_file_i=1:num_Nlg_files % for each Nlg .DAT file
            Nlg_file_number_string=file_start_details{Nlg_file_i}(end-2:end); % eg. find "003" from "File started. File index: 003"
            Nlg_file_name=[Nlg_file_name_letters '_' Nlg_file_number_string '.DAT']; % eg. "NEUR_003.DAT"
            if ~exist(fullfile(Nlg_folder,Nlg_file_name),'file')
                indices_missing_Nlg_files=[indices_missing_Nlg_files; Nlg_file_i]; %#ok<AGROW> % the indices here start from 1, unlike the Nlg. DAT file names which start from "000"
                disp(['Cannot open ' Nlg_file_name '; continuing to next file...']);
                continue
            end
            file_id=fopen(fullfile(Nlg_folder,Nlg_file_name)); % open the .DAT file for binary read access by MATLAB
            file_data=fread(file_id,Nlg_data_type); % import the Nlg AD count data as a single column vector with the default class "double"
            fclose(file_id); % close the .DAT file
            disp(['Reading Nlg file ' Nlg_file_number_string '...'])
            if first_file_loaded==1 % this is true for the first .DAT file loaded
                first_file_loaded=0; % this is set to zero after the first file has already been loaded
                if save_in_mat_or_Nlx_format==1 % if saving in MATLAB format
                    if isempty(samples_per_channel_per_file)
                        samples_per_channel_per_file=length(file_data)/num_channels;
                    end
                    AD_count_all_channels_all_files_int16=zeros(num_active_channels,samples_per_channel_per_file*num_Nlg_files,'int16'); % the variable where all voltage data will be saved, with the signed 16-bit integer format
                    timestamps_of_first_samples_usec=zeros(num_active_channels,num_Nlg_files); % the time stamps of the first sample of each file for each channel
                    % These are the only time stamps that will be saved if saving in
                    % the .mat format, because (1) these are the only time stamps the
                    % Nlg actually logs, (2) the time stamps of all samples can be
                    % calculated from these, and (3) this saves storage space. If we do
                    % want to save a time stamp for every sample, we should initialize
                    % here a "uint16" zeros vector the same size as
                    % "AD_count_all_channels_all_files_int16" above. The time stamps
                    % are microseconds from midnight, there are 8.64e10 microseconds in
                    % a day, so the unsigned 16-bit integer format is sufficient to
                    % store the time stamps, and saves storage and memory compared to
                    % other formats.
                    unwritten_samples=cell(size(partially_filled_Nlg_files)); % each cell contains the indices of the unwritten samples in a partially filled Nlg data file; the indices are counting from the beginning of the recording for a single channel and are the same for all channels
                    indices_of_first_samples=1:samples_per_channel_per_file:(num_Nlg_files-1)*samples_per_channel_per_file+1; % for a given channel, the indices of the first sample of every Nlg .DAT file, counting from the beginning of the recording; this is the same for all recording channels
                elseif save_in_mat_or_Nlx_format==2 % if saving in Nlx format
                    Nlx_512_block_sampling_period_usec=512*sampling_period_sec*1e6; % 512 times the sampling period, in us
                end
            end
            
            AD_count_data=reshape(file_data,num_channels,[]); % reshape the data: each row is the data for a channel
            AD_count_data=AD_count_data(active_channels,:); % take only the active channels
            if any(Nlg_file_i==partially_filled_Nlg_files) % if the current Nlg .DAT file is partially filled
                last_recorded_sample=find(diff(file_data),1,'last'); % the unwritten samples at the end of the file should all have differences between consecutive samples equal to 0, so the index of the last nonzero difference is the index of the last recorded sample
                if any(file_data(last_recorded_sample+1:end)~=Nlg_unwritten_data_value)
                    disp(['Error finding the unwritten samples in File ' Nlg_file_number_string ': the last consecutive run of data values does not equal the default value.'])
                end
                if mod(last_recorded_sample,num_channels)~=0 % checks for any remainder from the division of the index of the last recorded sample by the number of channels
                    disp(['Nlg did not finish recording from all channels in the last AD conversion period before stopping recording during File ' Nlg_file_number_string '; data from the unfinished AD conversion period discarded.'])
                    single_channel_first_unwritten_sample=floor(last_recorded_sample/num_channels)+1;
                else
                    single_channel_first_unwritten_sample=last_recorded_sample/num_channels+1;
                end
                if save_in_mat_or_Nlx_format==1 % if saving in MATLAB format
                    partially_filled_file_i=find(Nlg_file_i==partially_filled_Nlg_files);
                    unwritten_samples{partially_filled_file_i}=(single_channel_first_unwritten_sample:size(AD_count_data,2))+indices_of_first_samples(partially_filled_Nlg_files(partially_filled_file_i))-1;
                end
            end
            AD_count_data=AD_count_data-AD_count_for_zero_voltage; % subtract the AD count that represents 0 voltage
            
            if any(reference_channel) % if a channel will be used as reference
                non_reference_channels=active_channels~=reference_channel; % all the channels that are not the reference channel
                AD_count_data(non_reference_channels,:)=AD_count_data(non_reference_channels,:)-repmat(AD_count_data(active_channels==reference_channel,:),num_active_channels-1,1); % subtract the AD counts (equivalently, voltages) of the reference channel from those of all other channels
            end
            AD_count_data=-AD_count_data; % the voltage trace is inverted so that spikes start with an increase in voltage followed by a decrease
            
            file_start_timestamps_usec_from_sampling_period=file_start_timestamps_usec(first_file_after_recording_start)+samples_per_channel_per_file*(Nlg_file_i-first_file_after_recording_start)*sampling_period_sec*1e6; % time stamp of the first sample of the current .DAT file = the number of samples since recording start * sampling period
            file_start_timestamp_discrepancies_ms(Nlg_file_i)=round(file_start_timestamps_usec(Nlg_file_i)/1000)-round(file_start_timestamps_usec_from_sampling_period/1000); % difference between the event log time stamp and the time stamp calculated by counting samples, after roundig both to integer ms
            if abs(file_start_timestamp_discrepancies_ms(Nlg_file_i))>0
                disp(['File start time stamp read from event log and that calculated by samples differ by ' num2str(abs(file_start_timestamp_discrepancies_ms(Nlg_file_i))) ' ms...']) % occasional differences of 1 ms might be OK, because they could result from the millisecond resolution of the first time stamp after recording start
            end
            if timestamps_from_eventlog_or_sampling_period==1 || Nlg_file_i==first_file_after_recording_start % take the time stamp of the first sample of a .DAT file from the event log if the user specifies so, or if the current file is the first file after starting recording
                current_file_start_timestamps_usec=file_start_timestamps_usec(Nlg_file_i)+(0:ADC_sampling_period_usec:(num_channels-1)*ADC_sampling_period_usec); % time stamps of first sample of each of the channels; Jacob Vecht confirmed that the first sample of a file occurs at the file start time, unlike what was implied in the Neurologger manual (version 03-Apr-15 16:10:00)
            elseif timestamps_from_eventlog_or_sampling_period==2 % calculate the time stamp of the first sample of a .DAT file by counting the number of samples since recording start and multiplying by the sampling period
                current_file_start_timestamps_usec=file_start_timestamps_usec_from_sampling_period+(0:ADC_sampling_period_usec:(num_channels-1)*ADC_sampling_period_usec);
            end
            timestamps_of_first_samples_usec(:,Nlg_file_i)=current_file_start_timestamps_usec(active_channels); % only save the time stamps for the active channels
            
            if any(Nlg_file_i==partially_filled_Nlg_files) && Nlg_file_i~=num_Nlg_files % if the current Nlg .DAT file is partially filled and recording resumed later
                first_file_after_recording_start=Nlg_file_i+1; % updated to be the file after the partially filled file
            end
            
            if save_in_mat_or_Nlx_format==1 % if saving in MATLAB format
                AD_count_all_channels_all_files_int16(:,(Nlg_file_i-1)*samples_per_channel_per_file+1:Nlg_file_i*samples_per_channel_per_file)=int16(AD_count_data); % the data from the current .DAT file, converted to signed 16-bit integers
            elseif save_in_mat_or_Nlx_format==2 % if saving in Nlx format
                if any(Nlg_file_i==partially_filled_Nlg_files) % if the current .DAT file is partially filled
                    last_sample_in_channel_with_recorded_data=floor(last_recorded_sample/num_channels);
                else
                    last_sample_in_channel_with_recorded_data=samples_per_channel_per_file;
                end
                for active_channel_i=1:num_active_channels
                    AD_count_single_channel=AD_count_data(active_channel_i,1:last_sample_in_channel_with_recorded_data); % data from a single channel
                    AD_count_single_channel_blocks=vec2mat(AD_count_single_channel,512).'; % put the row vector of data into a 512-column matrix, and then transpose to get a 512-row matrix
                    num_blocks=size(AD_count_single_channel_blocks,2); % the number of 512-element blocks
                    timestamps_of_blocks_usec=current_file_start_timestamps_usec(active_channel_i)+(0:Nlx_512_block_sampling_period_usec:(num_blocks-1)*Nlx_512_block_sampling_period_usec); % the timestamps of the first sample in each 512-sample block, as required for exporting to the Nlx .ncs format
                    Nlx_file_name=['CSC' num2str(active_channels(active_channel_i)-1) '.ncs']; % Neuralynx numbers the first channel as channel 0, the second channel as channel 1, etc.
                    file_name_to_save=fullfile(output_folder,Nlx_file_name);
                    
                    % The following works with the current version of Mat2NlxEV
                    % (Version 6.0.0) as of this writing (7/6/2016).
                    AppendToFileFlag=1; % Append to the file if it already exists, otherwise create the file
                    ExportMode=1; % export all data
                    ExportModeVector=[]; % don't need this when exporting all data
                    FieldSelectionFlags=[1 0 0 0 1 0]; % whether to export: Timestamps, Channel Numbers, Sample Frequency, Number of Valid Samples, Samples, Header
                    Mat2NlxCSC(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,timestamps_of_blocks_usec,AD_count_single_channel_blocks);
                    % Note that the Nlx .ncs format only stores time in
                    % integer microseconds, and Mat2NlxCSC rounds
                    % non-integer times towards zero (instead of towards
                    % the nearest integer), so make sure all time stamps
                    % are integers before exporting (as they are in this
                    % code)
                end
            end
        end
        if save_in_mat_or_Nlx_format==1 && length(indices_missing_Nlg_files)~=num_Nlg_files % if saving in MATLAB format and at least one .DAT file was processed
            AD_count_all_channels_all_files_int16(:,[unwritten_samples{:}])=[]; % delete the unwritten samples
            for partially_filled_file_i=1:length(partially_filled_Nlg_files) % for each partially filled file, subtract the number of its unwritten samples from the indices of the first sample of all files after it
                files_after_partially_filled_file_i=partially_filled_Nlg_files(partially_filled_file_i)+1:num_Nlg_files;
                indices_of_first_samples(files_after_partially_filled_file_i)=indices_of_first_samples(files_after_partially_filled_file_i)-length(unwritten_samples{partially_filled_file_i});
            end
            
            clear variables_to_save
            for active_channel_i=1:num_active_channels
                file_name_to_save=fullfile(output_folder,['CSC' num2str(active_channels(active_channel_i)-1) '.mat']); % going with the Neuralynx numbering convention: the first channel is channel 0, the second channel is channel 1, etc.
                variables_to_save.AD_count_int16=AD_count_all_channels_all_files_int16(active_channel_i,:);
                variables_to_save.indices_of_first_samples=indices_of_first_samples;
                variables_to_save.timestamps_of_first_samples_usec=timestamps_of_first_samples_usec(active_channel_i,:);
                variables_to_save.samples_per_channel_per_file=samples_per_channel_per_file;
                variables_to_save.sampling_period_usec=sampling_period_sec*1e6;
                variables_to_save.AD_count_to_uV_factor=AD_count_to_uV_factor;
                save(file_name_to_save,'-struct','variables_to_save')
            end
        end
        disp(['Data from ' num2str(num_Nlg_files-length(indices_missing_Nlg_files)) ' out of ' num2str(num_Nlg_files) ' Nlg .DAT files in "' Nlg_folder '" were processed and saved.'])
    end
    %%
    % Save the options and parameters that were used
    if save_options_parameters_CD_figure
        file_name_to_save=fullfile(output_folder,['extract_Nlg_data_paramters_' date '.mat']);
        date_time_of_processing=datetime; % the date and time when this code was run
        clear variables_to_save
        variables_to_save.Nlg_folder=Nlg_folder;
        variables_to_save.output_folder=output_folder;
        variables_to_save.save_in_mat_or_Nlx_format=save_in_mat_or_Nlx_format;
        variables_to_save.save_event_file=save_event_file;
        variables_to_save.save_voltage_AD_count_files=save_voltage_AD_count_files;
        variables_to_save.inactive_channels=inactive_channels;
        variables_to_save.reference_channel=reference_channel;
        variables_to_save.global_or_local_clock_difference_estimation=global_or_local_clock_difference_estimation;
        variables_to_save.timestamps_from_eventlog_or_sampling_period=timestamps_from_eventlog_or_sampling_period;
        variables_to_save.Nlg_file_name_letters=Nlg_file_name_letters;
        variables_to_save.num_channels=num_channels;
        variables_to_save.AD_count_for_zero_voltage=AD_count_for_zero_voltage;
        variables_to_save.AD_count_to_uV_factor=AD_count_to_uV_factor;
        variables_to_save.sampling_period_sec=sampling_period_sec;
        variables_to_save.samples_per_channel_per_file=samples_per_channel_per_file;
        variables_to_save.Nlg_data_type=Nlg_data_type;
        variables_to_save.Nlg_unwritten_data_value=Nlg_unwritten_data_value;
        variables_to_save.old_clock_difference_bug_correction=old_clock_difference_bug_correction;
        variables_to_save.date_time_of_processing=date_time_of_processing;
        if exist('indices_missing_Nlg_files','var')
            variables_to_save.indices_missing_Nlg_files=indices_missing_Nlg_files;
        end
        if exist('file_start_timestamp_discrepancies_ms','var')
            variables_to_save.file_start_timestamp_discrepancies_ms=file_start_timestamp_discrepancies_ms;
        end
        variables_to_save.last_code_update=last_code_update;
        save(file_name_to_save,'-struct','variables_to_save')
    end
end
