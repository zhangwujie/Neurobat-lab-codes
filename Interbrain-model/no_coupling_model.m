% Comparing models with and without functional across-brain coupling, with
% the same behavioral statistics.
%%
% Parameters
num_simulations=100;

W_S=1; % functional self-coupling strength is -W_S
W_I=0.4; % functional across-brain coupling strength is W_I
tau_sec=15; % time constant in seconds
noise_type=1; % 0: no noise; 1: Gaussian white noise (uncorrelated over time and across bats; equal power across timescales) with zero mean and standard deviation "Gaussian_noise_std"
Gaussian_noise_std=0.15;
saved_behavior_transition_matrices='Behavior transition matrices.mat';
saved_inputs='Inputs, mean activity, 20200221.mat';
simulated_session_length_minute=100;
sampling_period=2.5012224;
%%
% Parameters for plotting
left_margin=0.09;
top_margin=0.1;
horizontal_margin=0.06;
vertical_margin=0.1;
activity_and_behavior_margin=0.006;
behavior_plot_height=0.045;
panel_width=0.4;
panel_height=0.3;
behavior_raster_height=0.9; % 1 is the maximum height
with_coupling_model_color=[162 87 255]/255;
without_coupling_model_color=[0 0.4470 0.7410];
bat_colors=[0 0.5 0; 1 0.6 0];
behavior_colors=...
    [51 51 204;... % Resting
    0 204 255;... % Active non-social
    200 142 255;... % Self grooming
    120 60 144;... % Social grooming
    0 255 0;... % Probing
    255 0 0;... % Fighting
    200 140 80;... % Mating
    153 153 0;... % Covering with wings
    234 234 0;... % Attempted reaching
    255 170 140;... % Blocking
    102 102 102;... % Other interactions
    0 0 0;... % Multiple behaviors
    ]/255;
font_size=6;
line_width=1;
activity_line_width=0.6;
legend_line_width=1.5;
tick_length=0.03;
y_range=[-0.5 1];
%%
num_time_points=round(simulated_session_length_minute*60/sampling_period);
sampling_period_minute=sampling_period/60;
interbrain_corr_model=nan(num_simulations,2,2); % simulation X full session/coordinated behavior removed X with/without across-bat coupling
%%
W_with_coupling=[-W_S W_I; W_I -W_S];
W_without_coupling=[-W_S 0; 0 -W_S];
%%
compare1=load(saved_behavior_transition_matrices,'behavior_group_names');
compare2=load(saved_inputs,'behavior_group_names');
if ~isequal(compare1.behavior_group_names,compare2.behavior_group_names)
    disp('Different behaviors for saved transition matrices and inputs.')
    return
end

load(saved_behavior_transition_matrices)
load(saved_inputs)

num_behaviors=length(behavior_group_names);
%%
all_simulated_activity=nan(2,num_time_points,num_simulations,2); % bat X time X simulation X with/without coupling
all_simulated_behaviors=nan(2,num_time_points,num_simulations); % bat X time X simulation; same behaviors are used for the with and without coupling simulations
noise_instantiations=nan(2,num_time_points,num_simulations); % bat X time X simulation; same inputs are used for the with and without coupling simulations
behavioral_statistics_from_one_or_two_chambers=1;
for simulation_i=1:num_simulations
    %%
    % Simulate behavior
    states_over_time=nan(num_time_points,1);
    states_over_time(1)=assign_state(initial_state_distribution{behavioral_statistics_from_one_or_two_chambers}); % use the initial distribution of the one-chamber sessions
    for step_i=2:num_time_points
        states_over_time(step_i)=assign_state(behavior_transition_matrix{behavioral_statistics_from_one_or_two_chambers}(:,states_over_time(step_i-1))); % use the transition matrix of the one-chamber sessions
    end
    %%
    simulated_behaviors=behaviors_for_each_state{behavioral_statistics_from_one_or_two_chambers}(:,states_over_time); % bat X time
    inputs=inputs_for_behaviors(simulated_behaviors);
    deterministic_input=inputs;
    if noise_type==1
        input_noise=randn(size(inputs))*Gaussian_noise_std;
        inputs=inputs+input_noise;
        noise_instantiations(:,:,simulation_i)=input_noise;
    end
    
    all_simulated_behaviors(:,:,simulation_i)=simulated_behaviors;
    
    times_with_coordinated_behavior=simulated_behaviors(1,:)==simulated_behaviors(2,:);
    %%
    for with_or_without_coupling=1:2
        if with_or_without_coupling==1
            W=W_with_coupling;
        elseif with_or_without_coupling==2
            W=W_without_coupling;
        end
        %%
        dynamics_equation=@(t,a) (W*a+interp1((1:num_time_points)*sampling_period,inputs',t,'linear','extrap')')/tau_sec;
        [~,simulated_activity]=ode45(dynamics_equation,(1:num_time_points)*sampling_period,-W\deterministic_input(:,1));
        
        all_simulated_activity(:,:,simulation_i,with_or_without_coupling)=simulated_activity';
        %%
        interbrain_corr_model(simulation_i,1,with_or_without_coupling)=corr(simulated_activity(:,1),simulated_activity(:,2));
        interbrain_corr_model(simulation_i,2,with_or_without_coupling)=corr(simulated_activity(~times_with_coordinated_behavior,1),simulated_activity(~times_with_coordinated_behavior,2));
    end
end
%%
figure
hold on
mean_to_plot=mean(interbrain_corr_model(:,:,1),1,'omitnan');
std_to_plot=std(interbrain_corr_model(:,:,1),[],1,'omitnan');
errorbar(1:2,mean_to_plot,std_to_plot,'-','Color',with_coupling_model_color)

mean_to_plot=mean(interbrain_corr_model(:,:,2),1,'omitnan');
std_to_plot=std(interbrain_corr_model(:,:,2),[],1,'omitnan');
errorbar(1:2,mean_to_plot,std_to_plot,'-','Color',without_coupling_model_color)

ylabel('Inter-brain correlation')
ylim(y_range)
xlim([0.5 2.5])
legend('Model: with functional across-brain coupling','Model: without functional across-brain coupling','location','eastoutside')
legend('boxoff')
set(gca,'XTick',1:2,'XTickLabel',{'Entire session' 'Coordinated behaviors removed'},'XTickLabelRotation',20,'TickLength',[tick_length 0.025],'box','off')
%%
example_simulation=1;
simulated_behaviors=all_simulated_behaviors(:,:,example_simulation);
times_with_coordinated_behavior=simulated_behaviors(1,:)==simulated_behaviors(2,:);
behavior_indices_over_time=simulated_behaviors';
behavior_uncoordinated=behavior_indices_over_time(~times_with_coordinated_behavior,:);
full_session_num_time_points=round(simulated_session_length_minute*60/sampling_period);
behavior_LBWH=cell(2,2); % with/without coupling X full session/coordinated behaviors removed
two_bat_activity_LBWH=cell(2,2); % with/without coupling X full session/coordinated behaviors removed
for with_or_without_coupling=1:2
    for full_or_part_session=1:2
        behavior_LBWH{with_or_without_coupling,full_or_part_session}=[left_margin+(panel_width+horizontal_margin)*(with_or_without_coupling-1) 1-top_margin-behavior_plot_height-(full_or_part_session-1)*(activity_and_behavior_margin+panel_height+vertical_margin+behavior_plot_height) panel_width behavior_plot_height];
        two_bat_activity_LBWH{with_or_without_coupling,full_or_part_session}=[behavior_LBWH{with_or_without_coupling,full_or_part_session}(1) behavior_LBWH{with_or_without_coupling,full_or_part_session}(2)-activity_and_behavior_margin-panel_height panel_width panel_height];
    end
end
figure('units','normalized','outerposition',[0 0 1 1]);
for with_or_without_coupling=1:2
    current_activity=squeeze(all_simulated_activity(:,:,example_simulation,:)); % bat X time X with/without coupling
    current_corr=corr(current_activity(1,:,with_or_without_coupling)',current_activity(2,:,with_or_without_coupling)');
    
    axes('Position',behavior_LBWH{with_or_without_coupling,1});
    plot_behavior_rasters(behavior_indices_over_time,sampling_period_minute,behavior_raster_height,behavior_colors)
    xlim([0 full_session_num_time_points-1]*sampling_period_minute)
    if with_or_without_coupling==1
        title({'Model: with functional across-brain coupling'; ['Correlation: ' num2str(round(current_corr*1000)/1000)]})
    elseif with_or_without_coupling==2
        title({'Model: without functional across-brain coupling'; ['Correlation: ' num2str(round(current_corr*1000)/1000)]})
    end
    set(gca,'XTick',[],'XColor','none','YColor','w','YTick',[behavior_raster_height/2 1+behavior_raster_height/2],'YTickLabel',{['\color[rgb]{' num2str(bat_colors(1,:)) '}Bat 1'] ['\color[rgb]{' num2str(bat_colors(2,:)) '}Bat 2']},'TickLength',[0 0])
    axis ij
    
    axes('Position',two_bat_activity_LBWH{with_or_without_coupling,1});
    hold on
    for trace_i=1:2
        plot((0:full_session_num_time_points-1)*sampling_period_minute,current_activity(trace_i,:,with_or_without_coupling),'Color',bat_colors(trace_i,:))
    end
    ylabel({'Neural activity' '(arbitrary unit)'})
    xlabel('Time (minute)')
    legend('Bat 1','Bat 2')
    xlim([0 full_session_num_time_points-1]*sampling_period_minute)
    
    % Coordinated behaviors removed
    current_corr=corr(current_activity(1,~times_with_coordinated_behavior,with_or_without_coupling)',current_activity(2,~times_with_coordinated_behavior,with_or_without_coupling)');
    
    axes('Position',behavior_LBWH{with_or_without_coupling,2});
    plot_behavior_rasters(behavior_uncoordinated,sampling_period_minute,behavior_raster_height,behavior_colors)
    xlim([0 size(behavior_uncoordinated,1)-1]*sampling_period_minute)
    title(['Correlation: ' num2str(round(current_corr*1000)/1000)])
    set(gca,'XTick',[],'XColor','none','YColor','w','YTick',[behavior_raster_height/2 1+behavior_raster_height/2],'YTickLabel',{['\color[rgb]{' num2str(bat_colors(1,:)) '}Bat 1'] ['\color[rgb]{' num2str(bat_colors(2,:)) '}Bat 2']},'TickLength',[0 0])
    axis ij
    
    axes('Position',two_bat_activity_LBWH{with_or_without_coupling,2});
    hold on
    for trace_i=1:2
        plot((0:size(behavior_uncoordinated,1)-1)*sampling_period_minute,current_activity(trace_i,~times_with_coordinated_behavior,with_or_without_coupling),'Color',bat_colors(trace_i,:))
    end
    ylabel({'Neural activity' '(arbitrary unit)'})
    xlabel('Time (minute)')
    legend('Bat 1','Bat 2')
    xlim([0 size(behavior_uncoordinated,1)-1]*sampling_period_minute)
end
%%
function state=assign_state(probabilities)
probability_cumulative_sums=cumsum(probabilities);
state_assignment_intervals=[[0; probability_cumulative_sums(1:end-1)] probability_cumulative_sums];
num_tries=0;
while true
    random_number=rand;
    state=find(random_number>state_assignment_intervals(:,1) & random_number<state_assignment_intervals(:,2));
    if length(state)==1
        break
    end
    num_tries=num_tries+1;
    if num_tries>20
        disp('Unable to assign state...')
        return
    end
end
end
%%
function plot_behavior_rasters(behavior_indices_over_time,sampling_period_minute,behavior_raster_height,behavior_colors)
hold on
for bat_i=1:2
    current_behavior=behavior_indices_over_time(1,bat_i);
    current_behavior_start_time=1;
    for time_i=1:size(behavior_indices_over_time,1)
        if behavior_indices_over_time(time_i,bat_i)~=current_behavior || time_i==size(behavior_indices_over_time,1)
            if any(current_behavior)
                rectangle_left=(current_behavior_start_time-1)*sampling_period_minute-sampling_period_minute/2;
                rectangle_bottom=bat_i-1;
                rectangle_width=((time_i-1)-current_behavior_start_time+1)*sampling_period_minute;
                rectangle('Position',[rectangle_left rectangle_bottom rectangle_width behavior_raster_height],'FaceColor',behavior_colors(current_behavior,:),'EdgeColor','none')
            end
            current_behavior=behavior_indices_over_time(time_i,bat_i);
            current_behavior_start_time=time_i;
        end
    end
end
ylim([0 1+behavior_raster_height])
end