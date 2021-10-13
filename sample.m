num_neuron = 100;
sample_rate = 100;%Hz
tot_time = 1;%s
num_time = tot_time*sample_rate;
baseline_FR = 50;
DOM_FR = 30;

rng('default');
PD = rand(num_neuron,1)*2*pi;
%PD = zeros(num_neuron,1);
x_(:,1)=cos(PD);
x_(:,2)=sin(PD);

movement_dir = [1,1];
movement_speed = gausswin(num_time);
movement_speed = movement_speed/sum(movement_speed);
movement_vel = movement_speed*movement_dir;
desired_FR = x_*movement_vel';
desired_FR = desired_FR./max(abs(desired_FR(:)));
desired_FR = baseline_FR + desired_FR*DOM_FR;

lamda = desired_FR/sample_rate;
spikes = poissrnd(lamda);

%% predict vel
norm_spikes = spikes - baseline_FR/sample_rate;
pred_Vel = x_'*norm_spikes;
pred_pos = cumsum(pred_Vel,2);

figure;plot(pred_pos','DisplayName','pred_pos');

figure;plot(pred_pos(1,:),pred_pos(2,:));axis equal

%% add noise/ reduce num_neurons and run through KF
% poisson noise
% white noise to lamda