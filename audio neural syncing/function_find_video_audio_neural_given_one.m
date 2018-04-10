% Given time points in one of the data streams (video, audio, or neural),
% find the corresponding time points in the other data streams. Need to
% have synced the data streams already. 5/17/2017 Wujie Zhang
%
% Inputs:
% -input_is_which_data: the string 'video', 'audio', or 'neural'
%
% -check_cameras_equal_num_frames: 0: don't check that the outputs for
% different cameras are consistent; any positive number: automatically
% adjust frame numbers if adjusting the frame numbers of one camera by less
% than this input would make the outputs for different cameras consistent
% (can enter a small number, eg. 0.01, to effectively ask every time
% whether to adjust or not)
%
% -input_times: video_frames_from_beginning, audio_times_ms_from_beginning,
% or neural_times_ms_from_midnight
%
% -neural_input_bat_name: if the input is neural times, specify here the
% name of the bat, can enter [] if neural data is only being recorded from
% one bat
%
% -video_input_camera_number: if the input is video frames, specify here
% the camera number, can enter [] if there is only one camera
%
% -video_neural_sync_method: for syncing video and neural data, 1 means
% fitting a single line over all points, 2 means interpolating between
% consecutive points; when there are very few skipped frames, use 1,
% otherwise use 2; defaults to 2 when it's empty
%
% -the other input variables are loaded from audio_neural_video_sync.mat
%
% Outputs:
% -neural_times_ms_from_midnight: the ith cell is the ith bat in
% "saved_bats_number_name"

function [video_frames_from_beginning,audio_times_ms_from_beginning,neural_times_ms_from_midnight]=function_find_video_audio_neural_given_one(input_is_which_data,check_cameras_equal_num_frames,input_times,neural_input_bat_name,video_input_camera_number,video_neural_sync_method,audio_sync_times_ms_from_beginning,main_neural_video_syncing_bat,neural_sync_with_audio_times_ms_from_midnight,neural_sync_with_video_times_ms_from_midnight,saved_bats_number_name,video_sync_frames_from_beginning,video_LED_on_off_frames,sync_events_indices)

if isequal(input_is_which_data,'video')
    video_frames_from_beginning=cell(length(video_sync_frames_from_beginning),1);
    video_frames_from_beginning{video_input_camera_number}=input_times;
    
    % Getting the neural times of the main video syncing bat
    neural_times_ms_from_midnight=cell(length(saved_bats_number_name),1);
    if length(saved_bats_number_name)==2
        index_main_neural_video_syncing_bat=find(ismember(saved_bats_number_name,main_neural_video_syncing_bat));
    elseif length(saved_bats_number_name)==1
        index_main_neural_video_syncing_bat=1;
    end
    if video_neural_sync_method==1
        [polyfit_neural_time_ms_from_midnight_as_function_of_frame.coefficients,~,polyfit_neural_time_ms_from_midnight_as_function_of_frame.mean_std]=polyfit(video_sync_frames_from_beginning{video_input_camera_number},neural_sync_with_video_times_ms_from_midnight{video_input_camera_number},1); % fitting the neural times of sync points as a linear function of the input frames
        neural_times_ms_from_midnight{index_main_neural_video_syncing_bat}=polyval(polyfit_neural_time_ms_from_midnight_as_function_of_frame.coefficients,video_frames_from_beginning{video_input_camera_number},[],polyfit_neural_time_ms_from_midnight_as_function_of_frame.mean_std);
    elseif isempty(video_neural_sync_method) || video_neural_sync_method==2
        neural_times_ms_from_midnight{index_main_neural_video_syncing_bat}=interp1(video_sync_frames_from_beginning{video_input_camera_number},neural_sync_with_video_times_ms_from_midnight{video_input_camera_number},video_frames_from_beginning{video_input_camera_number},'linear','extrap');
    end
    
    % Getting the audio times
    logger_index_convert_from=index_main_neural_video_syncing_bat;
    audio_times_ms_from_beginning=interp1(neural_sync_with_audio_times_ms_from_midnight{logger_index_convert_from},audio_sync_times_ms_from_beginning{logger_index_convert_from},neural_times_ms_from_midnight{logger_index_convert_from},'linear','extrap');
    
    % Getting the neural times from the other bat, when there is one
    if length(saved_bats_number_name)==2
        logger_index_convert_to=setdiff([1 2],logger_index_convert_from); % neural time stamps from this logger need to be calculated; this index indexes into the audio-neural syncing variables ("audio_times_ms_from_shared_pulse1", etc.)
        neural_times_ms_from_midnight{logger_index_convert_to}=interp1(audio_sync_times_ms_from_beginning{logger_index_convert_to},neural_sync_with_audio_times_ms_from_midnight{logger_index_convert_to},audio_times_ms_from_beginning,'linear','extrap'); % converting audio times to neural times of the other bat
    end
    
    % Getting video frames for other cameras
    for camera_i=1:length(video_sync_frames_from_beginning)
        if camera_i==video_input_camera_number
            continue
        end
        
        if video_neural_sync_method==1
            x_frames=reshape(video_LED_on_off_frames{video_input_camera_number}.',[],1);
            y_frames=reshape(video_LED_on_off_frames{camera_i}.',[],1);
            logical_indices_LED_transitions_not_shared=isnan(x_frames) | isnan(y_frames);
            x_frames(logical_indices_LED_transitions_not_shared)=[];
            y_frames(logical_indices_LED_transitions_not_shared)=[];
            [polyfit_frame_as_function_of_frame.coefficients,~,polyfit_frame_as_function_of_frame.mean_std]=polyfit(x_frames,y_frames,1); % fitting the frame number of LED transitions on the requested camera as a function of that on the input camera
            video_frames_from_beginning{camera_i}=polyval(polyfit_frame_as_function_of_frame.coefficients,video_frames_from_beginning{video_input_camera_number},[],polyfit_frame_as_function_of_frame.mean_std);
        elseif isempty(video_neural_sync_method) || video_neural_sync_method==2
            video_frames_from_beginning{camera_i}=interp1(neural_sync_with_video_times_ms_from_midnight{camera_i},video_sync_frames_from_beginning{camera_i},neural_times_ms_from_midnight{index_main_neural_video_syncing_bat},'linear','extrap');
        end
        if check_cameras_equal_num_frames>0
            if camera_i>1
                camera1_frame_minus_camera2_frame=mode(reshape(video_LED_on_off_frames{1}-video_LED_on_off_frames{2},[],1));
                if length(video_frames_from_beginning{1})>1 && ~isequal(diff(round(video_frames_from_beginning{camera_i})),diff(round(video_frames_from_beginning{1})))
                    if all((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames)
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: different numbers of frames requested for different cameras.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the numbers of frames requested for different cameras the same:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                elseif length(video_frames_from_beginning{1})==1 && round(round(video_frames_from_beginning{1})-round(video_frames_from_beginning{camera_i}))-camera1_frame_minus_camera2_frame~=0
                    if (video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: frame number difference between different cameras is not the same as the mode of the differences.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the frame number difference the same as the mode:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                end
            end
        end
    end
    
elseif isequal(input_is_which_data,'audio')
    audio_times_ms_from_beginning=input_times;
    
    % Getting the neural times
    neural_times_ms_from_midnight=cell(length(saved_bats_number_name),1);
    for bat_i=1:length(saved_bats_number_name)
        neural_times_ms_from_midnight{bat_i}=interp1(audio_sync_times_ms_from_beginning{bat_i},neural_sync_with_audio_times_ms_from_midnight{bat_i},audio_times_ms_from_beginning,'linear','extrap'); % converting audio times to neural times
    end
    
    % Getting the video frames
    video_frames_from_beginning=cell(length(video_sync_frames_from_beginning),1);
    if length(saved_bats_number_name)==2
        index_main_neural_video_syncing_bat=find(ismember(saved_bats_number_name,main_neural_video_syncing_bat));
    elseif length(saved_bats_number_name)==1
        index_main_neural_video_syncing_bat=1;
    end
    for camera_i=1:length(video_sync_frames_from_beginning)
        if video_neural_sync_method==1
            [polyfit_frame_as_function_of_neural_time_ms_from_midnight.coefficients,~,polyfit_frame_as_function_of_neural_time_ms_from_midnight.mean_std]=polyfit(neural_sync_with_video_times_ms_from_midnight{camera_i},video_sync_frames_from_beginning{camera_i},1); % fitting the frame numbers of sync points as a linear function of their neural times
            video_frames_from_beginning{camera_i}=polyval(polyfit_frame_as_function_of_neural_time_ms_from_midnight.coefficients,neural_times_ms_from_midnight{index_main_neural_video_syncing_bat},[],polyfit_frame_as_function_of_neural_time_ms_from_midnight.mean_std);
        elseif isempty(video_neural_sync_method) || video_neural_sync_method==2
            video_frames_from_beginning{camera_i}=interp1(neural_sync_with_video_times_ms_from_midnight{camera_i},video_sync_frames_from_beginning{camera_i},neural_times_ms_from_midnight{index_main_neural_video_syncing_bat},'linear','extrap');
        end
        if check_cameras_equal_num_frames>0
            if camera_i>1
                camera1_frame_minus_camera2_frame=mode(reshape(video_LED_on_off_frames{1}-video_LED_on_off_frames{2},[],1));
                if length(video_frames_from_beginning{1})>1 && ~isequal(diff(round(video_frames_from_beginning{camera_i})),diff(round(video_frames_from_beginning{1})))
                    if all((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames)
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: different numbers of frames requested for different cameras.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the numbers of frames requested for different cameras the same:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                elseif length(video_frames_from_beginning{1})==1 && round(round(video_frames_from_beginning{1})-round(video_frames_from_beginning{camera_i}))-camera1_frame_minus_camera2_frame~=0
                    if (video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: frame number difference between different cameras is not the same as the mode of the differences.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the frame number difference the same as the mode:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                end
            end
        end
    end
    
elseif isequal(input_is_which_data,'neural')
    neural_times_ms_from_midnight=cell(length(saved_bats_number_name),1);
    logger_index_convert_from=find(ismember(saved_bats_number_name,neural_input_bat_name)); % if neural data was recorded from two loggers, this logger's time stamps (the input) will be converted to those of the other logger; this index indexes into the audio-neural syncing variables ("audio_times_ms_from_shared_pulse1", etc.)
    neural_times_ms_from_midnight{logger_index_convert_from}=input_times;
    
    % Getting the audio times
    audio_times_ms_from_beginning=interp1(neural_sync_with_audio_times_ms_from_midnight{logger_index_convert_from},audio_sync_times_ms_from_beginning{logger_index_convert_from},neural_times_ms_from_midnight{logger_index_convert_from},'linear','extrap'); % converting neural times of the input bat to audio time
    
    % Getting the neural times from the other bat, when there is one
    if length(saved_bats_number_name)==2
        logger_index_convert_to=setdiff([1 2],logger_index_convert_from); % neural time stamps from this logger need to be calculated; this index indexes into the audio-neural syncing variables ("audio_times_ms_from_shared_pulse1", etc.)
        neural_times_ms_from_midnight{logger_index_convert_to}=interp1(audio_sync_times_ms_from_beginning{logger_index_convert_to},neural_sync_with_audio_times_ms_from_midnight{logger_index_convert_to},audio_times_ms_from_beginning,'linear','extrap'); % converting audio times to neural times of the other bat
    end
    
    % Getting video frames
    video_frames_from_beginning=cell(length(video_sync_frames_from_beginning),1);
    if length(saved_bats_number_name)==2
        index_main_neural_video_syncing_bat=find(ismember(saved_bats_number_name,main_neural_video_syncing_bat));
    elseif length(saved_bats_number_name)==1
        index_main_neural_video_syncing_bat=1;
    end
    for camera_i=1:length(video_sync_frames_from_beginning)
        if video_neural_sync_method==1
            [polyfit_frame_as_function_of_neural_time_ms_from_midnight.coefficients,~,polyfit_frame_as_function_of_neural_time_ms_from_midnight.mean_std]=polyfit(neural_sync_with_video_times_ms_from_midnight{camera_i},video_sync_frames_from_beginning{camera_i},1); % fitting the frame numbers of sync points as a linear function of their neural times
            video_frames_from_beginning{camera_i}=polyval(polyfit_frame_as_function_of_neural_time_ms_from_midnight.coefficients,neural_times_ms_from_midnight{index_main_neural_video_syncing_bat},[],polyfit_frame_as_function_of_neural_time_ms_from_midnight.mean_std);
        elseif isempty(video_neural_sync_method) || video_neural_sync_method==2
            video_frames_from_beginning{camera_i}=interp1(neural_sync_with_video_times_ms_from_midnight{camera_i},video_sync_frames_from_beginning{camera_i},neural_times_ms_from_midnight{index_main_neural_video_syncing_bat},'linear','extrap');
        end
        if check_cameras_equal_num_frames>0
            if camera_i>1
                camera1_frame_minus_camera2_frame=mode(reshape(video_LED_on_off_frames{1}-video_LED_on_off_frames{2},[],1));
                if length(video_frames_from_beginning{1})>1 && ~isequal(diff(round(video_frames_from_beginning{camera_i})),diff(round(video_frames_from_beginning{1})))
                    if all((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames)
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: different numbers of frames requested for different cameras.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the numbers of frames requested for different cameras the same:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                elseif length(video_frames_from_beginning{1})==1 && round(round(video_frames_from_beginning{1})-round(video_frames_from_beginning{camera_i}))-camera1_frame_minus_camera2_frame~=0
                    if (video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i}<check_cameras_equal_num_frames
                        video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                    else
                        adjust_or_not=questdlg({'Warning: frame number difference between different cameras is not the same as the mode of the differences.'; ['Frame adjustments for camera ' num2str(camera_i) ' that would make the frame number difference the same as the mode:']; num2str((video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame)-video_frames_from_beginning{camera_i})},'','Adjust','Don''t adjust','Adjust');
                        % the suggested adjustments are based on the lag between the two cameras; if the adjustment is one or a few frames, it should be ok, but if the adjustment is large, then it's probably wrong, maybe because one camera has skipped frames
                        if isequal(adjust_or_not,'Adjust')
                            video_frames_from_beginning{camera_i}=video_frames_from_beginning{1}-camera1_frame_minus_camera2_frame;
                            disp('Adjusted frame numbers.')
                        else
                            disp('Did not adjust frame numbers.')
                        end
                    end
                end
            end
        end
    end
end
