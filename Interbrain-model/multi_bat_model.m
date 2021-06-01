% Model of the neural activity of more than two interacting bats.
%%
% Parameters
num_simulations=100;
num_bats=4;
W_S=1; % functional self-coupling strength is -W_S
W_I=0.1; % functional across-brain coupling strength is W_I
tau_sec=15; % time constant in seconds

input_mean=0.2;
input_std=3.5;
input_moving_average_window_num_samples=1200;
noise_std=0.15;
num_diff_directions_to_sample=1000;
simulated_session_length_minute=100;
sampling_period=2.5;

plot_std_or_sem=1;
%%
S_weights=-W_S*eye(num_bats);

I_weights=W_I*ones(num_bats);
I_weights=I_weights-diag(diag(I_weights));

W=I_weights+S_weights;

sampling_freq=1/sampling_period;
num_time_points=round(simulated_session_length_minute*60/sampling_period);

moving_average_filter=ones(1,input_moving_average_window_num_samples)/input_moving_average_window_num_samples;
%%
mean_basis_vector=ones(num_bats,1)/num_bats;
basis_vector_length=norm(mean_basis_vector);
mean_unit_vector=mean_basis_vector/basis_vector_length;

I_matrix=eye(num_bats);
for col_i=1:num_bats
    I_matrix(:,col_i)=I_matrix(:,col_i)-(I_matrix(:,col_i)'*mean_unit_vector)*mean_unit_vector;
end
diff_unit_vectors=orth(I_matrix); % bat X basis vector; orthonormal basis vectors for the difference subspace
if size(diff_unit_vectors,2)~=num_bats-1
    disp('Error.')
    return
end
%%
logical_indices_across_brain_corrs=tril(true(num_bats),-1);
all_simulated_activity=nan(num_bats,num_time_points,num_simulations); % bat X time X simulation
all_inputs=nan(num_bats,num_time_points,num_simulations); % bat X time X simulation
power_spectral_centroid_mean_diff=nan(num_simulations,2); % simulation X mean/diff subspace
variance_mean_diff=nan(num_simulations,2); % simulation X mean/diff subspace
corr_bat_or_mean_diff=nan(num_simulations,2); % simulation X inter-brain corr / corr between mean and diff subspace

for simulation_i=1:num_simulations
    %%
    % Simulate behavior
    inputs=input_mean+randn(num_bats,num_time_points)*input_std;
    if any(input_moving_average_window_num_samples)
        for bat_i=1:num_bats
            inputs(bat_i,:)=cconv(inputs(bat_i,:),moving_average_filter,num_time_points);
        end
    end
    input_noise=randn(size(inputs))*noise_std;
    inputs=inputs+input_noise;
    %%
    dynamics_equation=@(t,a) (W*a+interp1((1:num_time_points)*sampling_period,inputs',t,'linear','extrap')')/tau_sec;
    [~,simulated_activity]=ode45(dynamics_equation,(1:num_time_points)*sampling_period,-W\inputs(:,1));
    %%
    mean_activity_projections=simulated_activity*mean_unit_vector;
    mean_activity=mean_activity_projections*mean_unit_vector'; % time X bat
    difference_subspace_activity=simulated_activity-mean_activity; % time X bat; activity in the mean direction projected out
    
    variance_mean_diff(simulation_i,1)=sum(var(mean_activity,1,1));
    variance_mean_diff(simulation_i,2)=sum(var(difference_subspace_activity,1,1))/(num_bats-1);
    
    current_signal=simulated_activity*mean_unit_vector;
    [power_spectrum,frequencies]=periodogram(current_signal-mean(current_signal),hamming(length(current_signal),'periodic'),[],sampling_freq,'onesided');
    power_spectral_centroid_mean_diff(simulation_i,1)=sum(power_spectrum.*frequencies)/sum(power_spectrum);
    
    Gaussian_random_numbers=randn(num_diff_directions_to_sample,num_bats-1);
    random_direction_PSCs=nan(num_diff_directions_to_sample,1);
    random_diff_mean_corr=nan(num_diff_directions_to_sample,1);
    for direction_i=1:num_diff_directions_to_sample
        random_direction=diff_unit_vectors*Gaussian_random_numbers(direction_i,:)';
        random_direction=random_direction/norm(random_direction);
        current_signal=simulated_activity*random_direction;
        [power_spectrum,frequencies]=periodogram(current_signal-mean(current_signal),hamming(length(current_signal),'periodic'),[],sampling_freq,'onesided');
        random_direction_PSCs(direction_i)=sum(power_spectrum.*frequencies)/sum(power_spectrum);
        
        random_diff_mean_corr(direction_i)=corr(current_signal,mean_activity_projections);
    end
    power_spectral_centroid_mean_diff(simulation_i,2)=mean(random_direction_PSCs);
    
    corr_matrix=corr(simulated_activity);
    corr_bat_or_mean_diff(simulation_i,1)=mean(corr_matrix(logical_indices_across_brain_corrs));
    corr_bat_or_mean_diff(simulation_i,2)=mean(random_diff_mean_corr);
    %%
    all_simulated_activity(:,:,simulation_i)=simulated_activity';
    all_inputs(:,:,simulation_i)=inputs;
end
%%
% Plot analysis results
ylabels={'Variance' 'Frequency (Hz)' 'Correlation'};
measure_titles={'Variance' 'Power spectral centroid' 'Correlation'};

measures_to_plot={variance_mean_diff power_spectral_centroid_mean_diff corr_bat_or_mean_diff};
figure
for measure_i=1:length(measures_to_plot)
    current_data=measures_to_plot{measure_i};
    means_to_plot=mean(current_data,1);
    current_std=std(current_data,[],1);
    current_sem=current_std/sqrt(num_simulations);
    if plot_std_or_sem==1
        error_to_plot=current_std;
    elseif plot_std_or_sem==2
        error_to_plot=current_sem;
    end
    
    subplot(2,3,measure_i)
    errorbar(1:2,means_to_plot,error_to_plot,'-ob')
    xlim([0.5 2+0.5])
    if measure_i==3
        set(gca,'XTick',1:2,'XTickLabel',{'Inter-brain' 'Between mean and diff. subspace'},'XTickLabelRotation',45)
    else
        set(gca,'XTick',1:2,'XTickLabel',{'Mean across all bats' 'Difference subspace'},'XTickLabelRotation',45)
    end
    ylabel(ylabels{measure_i})
    title(measure_titles{measure_i})
end