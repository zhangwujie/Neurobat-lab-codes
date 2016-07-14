% Given a folder containing Neurologger (Nlg) data from one recording
% session, extract the event and voltage data, and save in either the
% MATLAB .mat format, or in the Neuralynx (Nlx) .nev format (for the event
% file) and .ncs format (for the raw voltage traces represented by AD
% counts in bits) in another folder.
% -The Neurologger event file contains time stamps from the Neurologger and
% its transceiver, with clock drift between the two -- we correct the clock
% drift so that all time stamps are those of the Neurologger clock.
% -Can choose whether or not to save event and voltage data files
% -Before running, need to convert the Neurologger .NLE event file to .xlsx
% format.
% 7/2/2016, Wujie Zhang

% Beta testing version. Do not use this to seriously process data.
% -inputs and outputs 
% -make sure that all logger version specific parameters are in front
% instead of in the middle of the code
%%
% Options and parameters
Nlg_folder='C:\Users\phyllo\Documents\Maimon\ephys\yr2016_bat71319_robin_Nlg\neurologger_recording20160711\nlgformat\nlg1\';
output_folder='C:\Users\phyllo\Documents\Maimon\ephys\yr2016_bat71319_robin_Nlg\neurologger_recording20160711\nlxformat_wujie\';

save_in_mat_or_Nlx_format=1; % 1: save in .mat format; 2: save in .nev and .ncs formats
save_event_file=1; % 0: don't save any event file; 1: save one event file containing all events; 2: save one event file containing all events and one event file for each event type
save_voltage_AD_count_files=1;

inactive_channels=[]; % a vector containing numbers from 1 to the number of channels, indicating the disabled channels; enter an empty vector if all channels are active
reference_channel=[]; % specify a channel from 1 to num_channels as reference channel, whose AD counts (equivalently, voltages) will be subtracted from those of all other channels; leave empty to use ground instead of any recording channel as reference
Nlg_file_name_letters='NEUR';

% Parameters specific to the Neurologger used
num_channels=16; % number of recording channels
AD_count_for_zero_voltage=2048; % the AD count that represents zero volt, which will be subtracted from the recorded AD counts
sampling_period_sec=(512/15)/1e6; % specify here the exact sampling period in s for the signal in a given recording channel, or leave it empty to read the sampling period from the Nlg event file header; this is the number of channels times the sampling period of the AD converter
% For the Neurolog-16 that we use as of this writing (7/4/2016),
% 512/15 = 34.13333... us is exactly 16 times the sampling period of the AD
% converter. The sampling period from the header (34.133328 us) is slightly
% inaccurate (would result in a difference of ~1 ms over ~2 hours compared
% to using 512/15 us).
samples_per_channel_per_file=524288; % specify here the number of samples per channel per Nlg data file
% For the Neurolog-16 that we use as of this writing (7/6/2016), each file
% contain 2^24 bytes, and each sample takes 2 bytes, so the samples per
% channel per file is 2^24 bytes / (2 bytes * 16 channels) = 524288
Nlg_data_type='uint16'; % the data type of the AD counts in the Nlg .DAT files
% For the Neurolog-16 that we use as of this writing (7/11/2016), the data
% type is 16-bit unsigned integer.
Nlg_unwritten_data_value=65535; % specify here the default value of Nlg data samples when they haven't been written with actual recorded data
% When recording is stopped (at the end of a recording session or in the
% middle), the current .DAT file is partially filled with recorded data,
% and all samples in the unfilled portion of the file are equal to their
% default value, which is 65535 (the 16-bit unsigned integer with all bits
% being ones) in the version of Neurolog-16 used as of this writing
% (7/11/2016), and is 0 in some other versions of Neurologgers.

%%
% read Nlg event file
Nlg_event_file_xlsx_name=fullfile(Nlg_folder,'EVENTLOG.xlsx');
[~,~,imported_Nlg_event_xlsx]=xlsread(Nlg_event_file_xlsx_name); % imports the excel sheet as a cell array, importing numbers as numbers, words as strings, and empty cells as NaNs
string_or_number=zeros(size(imported_Nlg_event_xlsx,1),1);
for row_i=1:size(imported_Nlg_event_xlsx,1)
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
first_numeric_row=find(string_or_number==2,1);

continued_event_indices=find(string_or_number==1);
continued_event_indices=continued_event_indices(continued_event_indices>first_numeric_row);
for continued_i=continued_event_indices
    continued_from_i=find(string_or_number(1:continued_i)==2,1,'last'); % the last row before the current "...Continued" row that has a number, ie. the event from which the current "...Continued" row continues
    imported_Nlg_event_xlsx{continued_from_i,6}=[imported_Nlg_event_xlsx{continued_from_i,6} imported_Nlg_event_xlsx{continued_i,6}];
end
imported_Nlg_event_xlsx(continued_event_indices,:)=[];

if isempty(sampling_period_sec)
    Nlg_header_rows=1:first_numeric_row-2; % the row before the first numeric row is the column labels (eg. 'Time Stamp', 'Event Type', etc.), and the rows before that are the header
    Nlg_header=[imported_Nlg_event_xlsx{Nlg_header_rows,1}];
    string_before_Ts='ADC Period = ';
    string_after_Ts='us';
    Ts_start_index=strfind(Nlg_header,string_before_Ts)+length(string_before_Ts);
    Ts_end_index=strfind(Nlg_header,string_after_Ts)-1;
    sampling_period_sec=str2double(Nlg_header(Ts_start_index:Ts_end_index))/1e6;
    disp('Sampling period was not specified in this code; reading from Neurologger event file header instead...')
end

event_timestamps_usec=cell2mat(imported_Nlg_event_xlsx(first_numeric_row:end,3))*1000; % time stamps, converted from ms to us
event_timestamps_source=imported_Nlg_event_xlsx(first_numeric_row:end,4); % source of time stamps: logger or transceiver
event_types=imported_Nlg_event_xlsx(first_numeric_row:end,5);
event_types_and_details=strcat(imported_Nlg_event_xlsx(first_numeric_row:end,5),{'. '},imported_Nlg_event_xlsx(first_numeric_row:end,6)); % event type and details are concatenated together, with a period and space between them, eg. "File started. File index: 001"
%%
% Change all transceiver time stamps to the corresponding times of the
% logger clock (the transceiver clock and logger clock run at different
% speeds)
string_before_CD='CD=';
pre_CD_positions=strfind(event_types_and_details,string_before_CD);
PC_comments_indices=find(~cellfun(@isempty,pre_CD_positions));
clock_differences_sec=nan(length(PC_comments_indices),1);
logger_times=nan(length(PC_comments_indices),1);
for PC_comment_i=1:length(PC_comments_indices)
    current_row_i=PC_comments_indices(PC_comment_i);
    string_after_CD=event_types_and_details{current_row_i}(pre_CD_positions{current_row_i}+length(string_before_CD):end);
    position_CD_end=find(string_after_CD==' ',1)-1;
    clock_differences_sec(PC_comment_i)=str2double(string_after_CD(1:position_CD_end)); % clock differences are in s
    logger_times(PC_comment_i)=event_timestamps_usec(current_row_i);
end
%     clock_differences_sec=(clock_differences_sec*16)+0.03; % only use this line for data collected when transciever had error, ie. Mindy's data from bat #95, "Primus"
clock_differences_usec=clock_differences_sec*1e6; % convert from s to us

transceiver_times=logger_times-clock_differences_usec;
indices_events_with_transceiver_time=find(ismember(event_timestamps_source,'Transceiver')); % PC-generated comments (the points used to interpolate in-between) uses logger time
interpolated_clock_differences=interp1(transceiver_times,clock_differences_usec,event_timestamps_usec(indices_events_with_transceiver_time),'linear','extrap');
event_timestamps_usec(indices_events_with_transceiver_time)=event_timestamps_usec(indices_events_with_transceiver_time)+interpolated_clock_differences; % clock difference = logger - transceiver

event_timestamps_usec=round(event_timestamps_usec); % round all time stamps to integer microseconds

figure
hold on
plot((logger_times-event_timestamps_usec(1))/(1e6*60),clock_differences_usec/1e3,'ro')
plot((event_timestamps_usec(indices_events_with_transceiver_time)-event_timestamps_usec(1))/(1e6*60),interpolated_clock_differences/1e3,'b.')
legend('Recorded clock differences','Interpolated clock differences')
hold on
ylabel('Logger time - transceiver time (ms)')
xlabel('Logger time (minutes)')
%%
% save Nlx .nev format event file
if save_event_file
    mkdir(output_folder); % make the output folder
    if save_in_mat_or_Nlx_format==1
        file_name_to_save=fullfile(output_folder,'EVENTS.mat');
        save(file_name_to_save,'event_timestamps_usec','event_types_and_details')
        
        if save_event_file==2 % create sperate event file for each event category
            event_type_list=unique(event_types);
            for event_type_i=1:length(event_type_list)
                current_event_type=event_type_list{event_type_i};
                indices_of_events_of_this_type=ismember(event_types,current_event_type);
                file_name_to_save=fullfile(output_folder,['EVENTS_' current_event_type '.mat']);
                variables_to_save.event_timestamps_usec=event_timestamps_usec(indices_of_events_of_this_type);
                variables_to_save.event_types_and_details=event_types_and_details(indices_of_events_of_this_type);
                save(file_name_to_save,'-struct','variables_to_save')
            end
        end
    elseif save_in_mat_or_Nlx_format==2
        file_name_to_save=fullfile(output_folder,'EVENTS.nev');
        
        % The following works with the current version of Mat2NlxEV
        % (Version 6.0.0) as of this writing (7/5/2016).
        AppendToFileFlag=0; % save new file or overwrite existing file, but do not append to existing file
        ExportMode=1; % export all events
        ExportModeVector=[]; % don't need this when exporting all events
        FieldSelectionFlags=[1 0 0 0 1 0]; % whether to export: Timestamps, Event IDs, TTLs, Extras, Event Strings, Header; the Neurologger TTLs are treated as regular events, so don't export TTLs here
        Mat2NlxEV(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,event_timestamps_usec.',event_types_and_details);
        % For Mat2NlxEV 6.0.0, the "Timestamps" input needs to be 1 X N,
        % and the "Event Strings" input needs to be N X 1 (unlike what its
        % documentation specifies); note that the Nlx .nev format only
        % stores time in integer microseconds, and Mat2NlxEV rounds
        % non-integer times towards zero (instead of towards nearest
        % integer), so make sure all time stamps are integers before
        % exporting
        
        if save_event_file==2 % create sperate event file for each event category
            event_type_list=unique(event_types);
            for event_type_i=1:length(event_type_list)
                current_event_type=event_type_list{event_type_i};
                indices_of_events_of_this_type=ismember(event_types,current_event_type);
                file_name_to_save=fullfile(output_folder,['EVENTS_' current_event_type '.nev']);
                Mat2NlxEV(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,event_timestamps_usec(indices_of_events_of_this_type).',event_types_and_details(indices_of_events_of_this_type));
            end
        end
    end
end
%%
% Save the AD count values (in bits) encoding the raw voltage traces in the
% Neuralynx .ncs format
if save_voltage_AD_count_files
    mkdir(output_folder); % make the output folder
    didnt_read_all_Nlg_file=0;
    multi_channel_sampling_period_usec=(sampling_period_sec*1e6)/num_channels; % the Nlg AD converter samples each channel sequentially, so its sampling period is the sampling period of a single channel divided by the number of channels
    active_channels=1:num_channels;
    active_channels(inactive_channels)=[];
    num_active_channels=length(active_channels);
    
    indices_file_start_events=find(ismember(event_types,'File started')); % find the 'File started' events
    file_start_timestamps_usec=event_timestamps_usec(indices_file_start_events);
    file_start_details=event_types_and_details(indices_file_start_events);
    num_Nlg_files=length(file_start_timestamps_usec);
    
    indices_stop_recording_events=find(ismember(event_types_and_details,'Mode change. Stopped recording'));
    find_partially_filled_Nlg_files=zeros(length(indices_stop_recording_events),1);
    for stop_recording_i=1:length(indices_stop_recording_events)
        current_stop_recording_index=indices_stop_recording_events(stop_recording_i);
        find_partially_filled_Nlg_files(stop_recording_i)=max(indices_file_start_events(indices_file_start_events<current_stop_recording_index));
    end
    find_partially_filled_Nlg_files=unique(find_partially_filled_Nlg_files);
    found_event_strings=cell2mat(event_types_and_details(find_partially_filled_Nlg_files));
    partially_filled_Nlg_files=find(ismember(file_start_details,found_event_strings)); % the indices of the partially filled files; note these indices start from 1, unlike the actual file names of the Nlg .DAT files which start from 000

    if save_in_mat_or_Nlx_format==1
        AD_count_all_channels_all_files_int16=zeros(num_active_channels,samples_per_channel_per_file*num_Nlg_files,'int16');
        timestamps_of_first_samples_usec=zeros(num_active_channels,num_Nlg_files);
        % Note: If we want to save a time stamp for every sample, we should
        % initialize here a "uint16" zeros vector the same size as
        % "AD_count_data_int16" above. The time stamps are microseconds
        % from midnight, there are 8.64e10 microseconds in a day, so the
        % unsigned 16-bit integer format is sufficient to store the time
        % stamps, and saves storage and memory compared to the default
        % "double" class.
        unwritten_samples=cell(size(partially_filled_Nlg_files)); % each cell contains the indices of the unwritten samples in a Nlg data file, the indices are counting from the beginning of the recording for a single channel and are the same for all channels
        indices_of_first_samples=1:samples_per_channel_per_file:(num_Nlg_files-1)*samples_per_channel_per_file+1; % for a given channel, the indices of the first samples of every Nlg file, counting from the beginning of the recording; this is the same for all recording channels
    elseif save_in_mat_or_Nlx_format==2
        Nlx_512_block_sampling_period_usec=512*sampling_period_sec*1e6; % 512 times the sampling period, in us
    end
    
    for Nlg_file_i=1:num_Nlg_files % loop over the Nlg voltage data files
        Nlg_file_number_string=file_start_details{Nlg_file_i}(end-2:end);
        Nlg_file_name=[Nlg_file_name_letters '_' Nlg_file_number_string '.DAT'];
        if ~exist(fullfile(Nlg_folder,Nlg_file_name),'file')
            didnt_read_all_Nlg_file=1;
            disp(['Cannot open ' Nlg_file_name '; continuing to next file...']);
            continue
        end
        file_id=fopen(fullfile(Nlg_folder,Nlg_file_name));
        file_data=fread(file_id,Nlg_data_type); % imports the Nlg data as a column vector with the default class "double"
        fclose(file_id);
        disp(['Reading Nlg file ' Nlg_file_number_string '...'])
        AD_count_data=reshape(file_data,num_channels,[]); % reshape the data: each row is the data for a channel
        AD_count_data=AD_count_data(active_channels,:);
        if any(Nlg_file_i==partially_filled_Nlg_files)
            if save_in_mat_or_Nlx_format==1
                partially_filled_file_i=find(Nlg_file_i==partially_filled_Nlg_files);
                unwritten_samples{partially_filled_file_i}=find(AD_count_data(1,:)==Nlg_unwritten_data_value)+indices_of_first_samples(partially_filled_Nlg_files(partially_filled_file_i))-1; % only checks one channel, because we assume that Nlg always finishes recording from all channels in the last AD conversion period before stopping recording, which seems to be the case (Jacob Vecht is unsure either; the next lines check for this)
            end
            last_sample_with_recorded_data=find(file_data==Nlg_unwritten_data_value,1)-1;
            if mod(last_sample_with_recorded_data,num_channels)~=0 % checks for any remainder from the division of the index of the last sample with actual recorded data, by the number of channels
                disp(['Nlg did not finish recording from all channels in the last AD conversion period before stopping recording during File ' Nlg_file_number_string '.'])
            end
        end
        AD_count_data=AD_count_data-AD_count_for_zero_voltage; % subtracts the AD count that represents 0 voltage
        
        if any(reference_channel)
            non_reference_channels=active_channels~=reference_channel; % all the channels that are not the reference channel
            AD_count_data(non_reference_channels,:)=AD_count_data(non_reference_channels,:)-repmat(AD_count_data(active_channels==reference_channel,:),num_active_channels-1,1); % subtract the AD counts (equivalently, voltages) of the reference channel from those of all other channels
        end
        AD_count_data=-AD_count_data; % the voltage trace is inverted to conform to convention
        
        current_file_start_timestamps_usec=file_start_timestamps_usec(Nlg_file_i)+(0:multi_channel_sampling_period_usec:(num_channels-1)*multi_channel_sampling_period_usec); % time stamps of first sample of each of the channels; Jacob Vecht confirmed that the first sample of a file occurs at the file start time, unlike what was implied in the Neurologger manual (version 03-Apr-15 16:10:00)
        current_file_start_timestamps_usec=current_file_start_timestamps_usec(active_channels);
        timestamps_of_first_samples_usec(:,Nlg_file_i)=current_file_start_timestamps_usec;
        if save_in_mat_or_Nlx_format==1
            AD_count_all_channels_all_files_int16(:,(Nlg_file_i-1)*samples_per_channel_per_file+1:Nlg_file_i*samples_per_channel_per_file)=int16(AD_count_data);
        elseif save_in_mat_or_Nlx_format==2
            if any(Nlg_file_i==partially_filled_Nlg_files)
                last_sample_in_channel_with_recorded_data=floor(last_sample_with_recorded_data/num_channels);
            else
                last_sample_in_channel_with_recorded_data=samples_per_channel_per_file;
            end
            for active_channel_i=1:num_active_channels % loop over the active recording channels
                AD_count_single_channel=AD_count_data(active_channel_i,1:last_sample_in_channel_with_recorded_data); % data from a single channel
                AD_count_single_channel_blocks=vec2mat(AD_count_single_channel,512).'; % put the row vector of data into a 512-column matrix, and then transpose to get a 512-row matrix
                num_blocks=size(AD_count_single_channel_blocks,2); % the number of 512-element blocks
                timestamps_of_blocks_usec=current_file_start_timestamps_usec(active_channel_i)+(0:Nlx_512_block_sampling_period_usec:(num_blocks-1)*Nlx_512_block_sampling_period_usec); % the timestamps of the first sample in each 512-sample block, as required for exporting to the Nlx .ncs format
                Nlx_file_name=['CSC' num2str(active_channels(active_channel_i)-1) '.ncs']; % Neuralynx numbers the first channel as channel 0, the second channel as channel 1, etc.
                file_name_to_save=fullfile(output_folder,Nlx_file_name);
                
                % The following works with the current version of Mat2NlxEV
                % (Version 6.0.0) as of this writing (7/6/2016).
                AppendToFileFlag=1;
                ExportMode=1;
                ExportModeVector=[];
                FieldSelectionFlags=[1 0 0 0 1 0];
                Mat2NlxCSC(file_name_to_save,AppendToFileFlag,ExportMode,ExportModeVector,FieldSelectionFlags,timestamps_of_blocks_usec,AD_count_single_channel_blocks);
                % Note that the Nlx .ncs format only stores time in integer
                % microseconds, and Mat2NlxCSC rounds non-integer times towards
                % zero (instead of towards nearest integer), so make sure all
                % time stamps are integers before exporting (as they are in
                % this code)
            end
        end
    end
    if save_in_mat_or_Nlx_format==1
        AD_count_all_channels_all_files_int16(:,[unwritten_samples{:}])=[];
        for partially_filled_file_i=1:length(partially_filled_Nlg_files)
            files_after_partially_filled_file_i=partially_filled_Nlg_files(partially_filled_file_i)+1:num_Nlg_files;
            indices_of_first_samples(files_after_partially_filled_file_i)=indices_of_first_samples(files_after_partially_filled_file_i)-length(unwritten_samples{partially_filled_file_i});
        end

        clear variables_to_save
        for active_channel_i=1:num_active_channels % loop over the active recording channels
            file_name_to_save=fullfile(output_folder,['CSC' num2str(active_channels(active_channel_i)-1) '.mat']); % going with the Neuralynx numbering convention: the first channel is channel 0, the second channel is channel 1, etc.
            variables_to_save.AD_count_int16=AD_count_all_channels_all_files_int16(active_channel_i,:);
            variables_to_save.indices_of_first_samples=indices_of_first_samples;
            variables_to_save.timestamps_of_first_samples_usec=timestamps_of_first_samples_usec(active_channel_i,:);
            variables_to_save.samples_per_channel_per_file=samples_per_channel_per_file;
            variables_to_save.sampling_period_usec=sampling_period_sec*1e6;
            save(file_name_to_save,'-struct','variables_to_save')
        end
    end
    if didnt_read_all_Nlg_file
        disp('Finished saving data from some but not all Nlg .DAT files.')
    else
        disp('Finished saving data from all Nlg .DAT files!')
    end
end