% Kuramoto model of the neural activity of two interacting bats.
%%
% Parameters
num_simulations_one_two_chambers=[100 100]; % [one-chamber two-chambers]

timescale_parameter=0.05;
natural_frequencies_mean=timescale_parameter*[0.1 0.2]; % the two natural frequencies are drawn lognormal distributions with these means
natural_frequencies_std=timescale_parameter*0.01; % the two natural frequencies are drawn lognormal distributions with this standard deviation
coupling_one_two_chambers=timescale_parameter*[0.07 0];
simulated_session_length_minute=100;
sampling_period=2.5012224;

plot_example_simulations=1;
%%
% Parameters for plotting
bat_colors=[0 0.5 0; 1 0.6 0];
activity_plot_left=0.04;
activity_plot_width=0.83;
activity_plot_height=0.31;
activity_plot_bottom=0.51;
vertical_margin=0.1;
font_size=15;
line_width=1.5;
bar_plot_line_width=1.5;
simulation_summary_color=[162 87 255]/255;
bar_plot_y_ranges=[0 10; 0 1.7]; % magnitude/timescale X y lower/upper limit
tick_length=0.03;
%%
mu=log((natural_frequencies_mean.^2)./sqrt(natural_frequencies_std^2+natural_frequencies_mean.^2));
sigma=sqrt(log(natural_frequencies_std^2./(natural_frequencies_mean.^2)+1));
%%
sampling_freq=1/sampling_period;
num_time_points=round(simulated_session_length_minute*60/sampling_period);
%%
all_simulated_activity=cell(2,1); % one/two-chambers
simulation_variance=cell(2,1); % one/two-chambers
simulation_power_spectral_centroid=cell(2,1); % one/two-chambers
natural_frequencies=cell(2,1); % one/two-chambers
for one_or_two_chamber=1:2
    all_simulated_activity{one_or_two_chamber}=nan(2,num_time_points,num_simulations_one_two_chambers(one_or_two_chamber)); % bat X time X simulation
    simulation_power_spectral_centroid{one_or_two_chamber}=nan(num_simulations_one_two_chambers(one_or_two_chamber),2); % simulation X mean/diff
    simulation_variance{one_or_two_chamber}=nan(num_simulations_one_two_chambers(one_or_two_chamber),2); % simulation X mean/diff
    natural_frequencies{one_or_two_chamber}=nan(2,num_simulations_one_two_chambers(one_or_two_chamber)); % bat 1/bat 2 X simulation
    for simulation_i=1:num_simulations_one_two_chambers(one_or_two_chamber)
        %%
        for bat_i=1:2
            natural_frequencies{one_or_two_chamber}(bat_i,simulation_i)=lognrnd(mu(bat_i),sigma(bat_i));
        end
        %%
        simulated_activity=simulate_modified_Kuramoto_model((1:num_time_points)*sampling_period,[0 0]',natural_frequencies{one_or_two_chamber}(:,simulation_i),coupling_one_two_chambers(one_or_two_chamber)); % initial condition is [0 0]'
        %%
        % Analyze the mean and difference
        mean_diff_traces=[(simulated_activity(:,1)+simulated_activity(:,2))/2 (simulated_activity(:,1)-simulated_activity(:,2))/2]; % time X mean/differenece
        %%
        for trace_i=1:2
            current_signal=mean_diff_traces(:,trace_i);
            
            [power_spectrum,frequencies]=periodogram(current_signal-mean(current_signal),hamming(length(current_signal),'periodic'),[],sampling_freq,'onesided');
            simulation_power_spectral_centroid{one_or_two_chamber}(simulation_i,trace_i)=sum(power_spectrum.*frequencies)/sum(power_spectrum);
            simulation_variance{one_or_two_chamber}(simulation_i,trace_i)=var(current_signal,1);
        end
        %%
        all_simulated_activity{one_or_two_chamber}(:,:,simulation_i)=simulated_activity';
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
        legend('Bat 1','Bat 2')
        xlabel('Time (minute)')
        ylabel('Activity')
        title(one_two_chamber_text{one_or_two_chamber})
        set(gca,'FontSize',font_size,'FontWeight','bold','LineWidth',line_width)
        
        subplot('Position',[activity_plot_left activity_plot_bottom-vertical_margin-activity_plot_height activity_plot_width activity_plot_height]);
        simulated_activity=all_simulated_activity{one_or_two_chamber}(:,:,simulation_i)';
        mean_difference=[(simulated_activity(:,1)+simulated_activity(:,2))/2 (simulated_activity(:,1)-simulated_activity(:,2))/2];
        hold on
        for bat_i=1:2
            plot(time_minutes,mean_difference(:,bat_i),'LineWidth',line_width)
        end
        legend('Mean','Difference')
        xlim(time_minutes([1 end]))
        xlabel('Time (minute)')
        ylabel('Activity')
        set(gca,'FontSize',font_size,'FontWeight','bold','LineWidth',line_width)
    end
end
%%
simulation_magnitude_timescale_measures={simulation_variance simulation_power_spectral_centroid};
measure_titles={{'Variance ratio,'; 'mean/difference'} {'Power spectral centroid ratio,'; 'mean/difference'}};

figure
for magnitude_or_timescale=1:2
    simulation_ratio_means=nan(2,1);
    simulation_ratio_stds=nan(2,1);
    simulation_current_ratios=cell(2,1);
    for one_two_chamber=1:2
        simulation_current_ratios{one_two_chamber}=simulation_magnitude_timescale_measures{magnitude_or_timescale}{one_two_chamber}(:,1)./simulation_magnitude_timescale_measures{magnitude_or_timescale}{one_two_chamber}(:,2);
        simulation_ratio_means(one_two_chamber)=mean(simulation_current_ratios{one_two_chamber});
        simulation_ratio_stds(one_two_chamber)=std(simulation_current_ratios{one_two_chamber},1);
    end
    
    subplot(1,2,magnitude_or_timescale)
    hold on
    errorbar(1:2,simulation_ratio_means,simulation_ratio_stds,'-','color',simulation_summary_color,'LineWidth',bar_plot_line_width)
    xlim([0.5 2.5])
    if any(bar_plot_y_ranges(magnitude_or_timescale,:))
        ylim(bar_plot_y_ranges(magnitude_or_timescale,:))
    end
    ylabel('Mean/difference ratio')
    title(measure_titles{magnitude_or_timescale})
    set(gca,'XTick',1:2,'XTickLabel',{'One-chamber' 'Two-chambers'},'XTickLabelRotation',20,'LineWidth',line_width,'FontSize',font_size,'FontWeight','bold','TickLength',[tick_length 0.025],'box','off')
end
%%
function simulated_activity=simulate_modified_Kuramoto_model(time_vector,phase_initial_condition,natural_frequencies,coupling)
num_bats=length(phase_initial_condition);
[~,simulated_phase]=ode15s(@dynamics_equation,time_vector,phase_initial_condition); % ode45 seems not to work
simulated_activity=(sin(simulated_phase)+1)/2;

    function dtheta_dt=dynamics_equation(t,theta)
        dtheta_dt=natural_frequencies+coupling*sum(sin(repmat(theta,1,num_bats)'-repmat(theta,1,num_bats)),2);
    end
end