% Main model of the neural activity of two interacting bats.
%%
% Parameters
num_simulations_one_two_chambers=[100 100]; % [one-chamber two-chambers]
W_S=1; % functional self-coupling strength is -W_S
W_I=0.4; % functional across-brain coupling strength is W_I
tau_sec=15; % time constant in seconds
noise_type=1; % 0: no noise; 1: Gaussian white noise with zero mean and standard deviation "Gaussian_noise_std"
Gaussian_noise_std=0.15;
saved_behavior_transition_matrices='Behavior transition matrices.mat';
saved_inputs='Inputs, mean activity, 20200221.mat';
simulated_session_length_minute=100;
sampling_period=2.5012224;

plot_example_simulations=1;
plot_std_or_sem=1;
%%
% Parameters for plotting
bat_colors=[0 0.5 0; 1 0.6 0];
behavior_raster_height=0.9; % 1 is the maximum height
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
    ]/255;
activity_plot_left=0.04;
activity_plot_width=0.83;
activity_plot_height=0.31;
activity_plot_bottom=0.51;
behavior_activity_margin=0.01;
behavior_plot_height=0.11;
font_size=15;
line_width=1.5;
%%
W_one_chamber=[-W_S W_I; W_I -W_S];
W_two_chambers=[-W_S 0; 0 -W_S];

sampling_freq=1/sampling_period;
num_time_points=round(simulated_session_length_minute*60/sampling_period);
%%
load(saved_behavior_transition_matrices)
load(saved_inputs)

compare1=load(saved_behavior_transition_matrices,'behavior_group_names');
compare2=load(saved_inputs,'behavior_group_names');
if ~isequal(compare1.behavior_group_names,compare2.behavior_group_names)
    disp('Different behaviors for saved transition matrices and inputs.')
    return
end
num_behaviors=length(behavior_group_names);
%%
all_simulated_activity=cell(2,1); % one/two-chambers
all_simulated_behaviors=cell(2,1); % one/two-chambers
power_spectral_centroid=cell(2,1); % one/two-chambers
variances=cell(2,1); % one/two-chambers
for one_or_two_chamber=1:2
    all_simulated_activity{one_or_two_chamber}=nan(2,num_time_points,num_simulations_one_two_chambers(one_or_two_chamber)); % bat X time X simulation
    all_simulated_behaviors{one_or_two_chamber}=nan(2,num_time_points,num_simulations_one_two_chambers(one_or_two_chamber)); % bat X time X simulation
    power_spectral_centroid{one_or_two_chamber}=nan(num_simulations_one_two_chambers(one_or_two_chamber),2); % simulation X mean/diff
    variances{one_or_two_chamber}=nan(num_simulations_one_two_chambers(one_or_two_chamber),2); % simulation X mean/diff
    if one_or_two_chamber==1
        W=W_one_chamber;
    elseif one_or_two_chamber==2
        W=W_two_chambers;
    end
    for simulation_i=1:num_simulations_one_two_chambers(one_or_two_chamber)
        %%
        % Simulate behavior
        states_over_time=nan(num_time_points,1);
        states_over_time(1)=assign_state(initial_state_distribution{one_or_two_chamber});
        for step_i=2:num_time_points
            states_over_time(step_i)=assign_state(behavior_transition_matrix{one_or_two_chamber}(:,states_over_time(step_i-1)));
        end
        %%
        simulated_behaviors=behaviors_for_each_state{one_or_two_chamber}(:,states_over_time); % bat X time
        inputs=inputs_for_behaviors(simulated_behaviors);
        deterministic_input=inputs;
        if noise_type==1
            input_noise=randn(size(inputs))*Gaussian_noise_std;
            inputs=inputs+input_noise;
        end
        %%
        dynamics_equation=@(t,a) (W*a+interp1((1:num_time_points)*sampling_period,inputs',t,'linear','extrap')')/tau_sec;
        [~,simulated_activity]=ode45(dynamics_equation,(1:num_time_points)*sampling_period,-W\deterministic_input(:,1));
        %%
        % Analyze the mean and difference
        activity_traces=[(simulated_activity(:,1)+simulated_activity(:,2))/2 (simulated_activity(:,1)-simulated_activity(:,2))/2]; % time X mean/differenece
       
        for trace_i=1:2
            current_signal=activity_traces(:,trace_i);
            
            [power_spectrum,frequencies]=periodogram(current_signal-mean(current_signal),hamming(length(current_signal),'periodic'),[],sampling_freq,'onesided');
            power_spectral_centroid{one_or_two_chamber}(simulation_i,trace_i)=sum(power_spectrum.*frequencies)/sum(power_spectrum);
            
            variances{one_or_two_chamber}(simulation_i,trace_i)=var(current_signal,1);
        end
        %%
        all_simulated_activity{one_or_two_chamber}(:,:,simulation_i)=simulated_activity';
        all_simulated_behaviors{one_or_two_chamber}(:,:,simulation_i)=simulated_behaviors;
    end
end
%%
% Plot example session
if plot_example_simulations
    one_two_chamber_text={'Example one-chamber simulation' 'Example two-chambers simulation'};
    time_minutes=(1:num_time_points)*sampling_period/60;
    sampling_period_minute=sampling_period/60;
    for one_or_two_chamber=1:2
        simulation_i=randi(num_simulations_one_two_chambers(one_or_two_chamber));
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot('Position',[activity_plot_left activity_plot_bottom activity_plot_width activity_plot_height]);
        
        hold on
        for bat_i=1:2
            plot(time_minutes,all_simulated_activity{one_or_two_chamber}(bat_i,:,simulation_i),'Color',bat_colors(bat_i,:),'LineWidth',line_width)
        end
        xlim(time_minutes([1 end]))
        xlabel('Time (minute)')
        ylabel('Activity')
        set(gca,'FontSize',font_size,'FontWeight','bold','LineWidth',line_width)
        
        behavior_plot_handle=subplot('Position',[activity_plot_left activity_plot_bottom+activity_plot_height+behavior_activity_margin activity_plot_width behavior_plot_height]);
        hold on
        for step_i=1:num_time_points
            for bat_i=1:2
                rectangle_left=step_i*sampling_period_minute-sampling_period_minute/2;
                rectangle_bottom=bat_i-1;
                rectangle('Position',[rectangle_left rectangle_bottom sampling_period_minute behavior_raster_height],'FaceColor',behavior_colors(all_simulated_behaviors{one_or_two_chamber}(bat_i,step_i,simulation_i),:),'EdgeColor','none')
            end
        end
        ylim([0 1+behavior_raster_height])
        xlim([1 num_time_points]*sampling_period_minute)
        set(gca,'XTick',[],'XColor','none','YColor','w','YTick',[behavior_raster_height/2 1+behavior_raster_height/2],'YTickLabel',{['\color[rgb]{' num2str(bat_colors(1,:)) '}Bat 1'] ['\color[rgb]{' num2str(bat_colors(2,:)) '}Bat 2']},'TickLength',[0 0],'FontSize',font_size,'FontWeight','bold','LineWidth',line_width)
        axis ij
        title(one_two_chamber_text{one_or_two_chamber})
        
        subplot('Position',[activity_plot_left+activity_plot_width behavior_plot_handle.Position(2)-activity_plot_height 1-(activity_plot_left+activity_plot_width) activity_plot_height+behavior_plot_handle.Position(4)])
        for category_i=1:length(behavior_group_names)
            rectangle('Position',[0 category_i 1 0.7],'FaceColor',behavior_colors(category_i,:),'EdgeColor','none')
            text(1.2,category_i+0.3,behavior_group_names{category_i},'FontSize',font_size*1.1,'FontWeight','bold')
        end
        axis([-0.3 10 1 category_i+0.7])
        axis ij
        axis off
    end
end
%%
% Plot analysis results
measures_to_plot={power_spectral_centroid variances};
measure_titles={{'Power spectral centroid ratio,'; 'mean/difference'} {'Variance ratio,'; 'mean/difference'}};

figure
for measure_i=1:length(measures_to_plot)
    ratio_means=nan(2,1); % one/two-chamber
    ratio_stds=nan(2,1);
    ratio_sems=nan(2,1);
    for one_or_two_chamber=1:2
        current_data=measures_to_plot{measure_i}{one_or_two_chamber};
        current_ratios=current_data(:,1)./current_data(:,2);
        
        ratio_means(one_or_two_chamber)=mean(current_ratios);
        ratio_stds(one_or_two_chamber)=std(current_ratios);
        ratio_sems(one_or_two_chamber)=std(current_ratios)/sqrt(length(current_ratios));
    end
    
    if plot_std_or_sem==1
        error_to_plot=ratio_stds;
    elseif plot_std_or_sem==2
        error_to_plot=ratio_sems;
    end
    
    subplot(2,2,3-measure_i)
    hold on
    errorbar(1:2,ratio_means,error_to_plot,'-o')
    xlim([0.5 2.5])
    set(gca,'XTick',1:2,'XTickLabel',{'One-chamber' 'Two-chambers'})
    ylabel('Mean/difference ratio')
    title(measure_titles{measure_i})
end
%%
% Simulating behavior using a Markov chain
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