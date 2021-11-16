n_neurons = 50;
n_tr = 2000;
dom = 5;
baseline_firing = 20;
neural_noise = 2;
init_pd = rand(n_neurons,1)*360;
init_dom = dom*rand(n_neurons,1);

init_tuning = init_dom.*[cosd(init_pd), sind(init_pd)];

init_reach_angle =  rand(n_tr,1)*360; 
init_reach_mag = 1*ones(n_tr,1);
init_reach_xy = init_reach_mag.*[cosd(init_reach_angle), sind(init_reach_angle)];

init_fr = baseline_firing + init_reach_xy*init_tuning';
init_fr = init_fr + neural_noise*randn(size(init_fr));

figure
scatter(init_reach_angle,init_fr(:,2))


[coeff,score,latent,~,~,mu] = pca(init_fr);

figure; plot(latent)
corr_mag = 0.5;

corr_reach_angle =  rand(n_tr,1)*360;
corr_reach_mag = corr_mag*ones(n_tr,1);
corr_reach_xy = corr_reach_mag.*[cosd(corr_reach_angle), sind(corr_reach_angle)];

%%
%Change how corrective tuning is related to initial
% corr_tuning = init_tuning;
%if corr Dom change
% corr_tuning = 0.5*init_tuning;
%if corr have neuron change
corr_neurons = 0.5*n_neurons;
corr_pd = rand(corr_neurons,1)*360;
corr_dom = dom*rand(corr_neurons,1);
corr_tuning = corr_dom.*[cosd(corr_pd), sind(corr_pd)];

corr_fr = baseline_firing + corr_reach_xy*corr_tuning';
corr_fr = corr_fr + neural_noise*randn(size(corr_fr));

corr_score = (corr_fr-mu)*coeff;


figure
scatter(score(:,1), score(:,2))
hold on
scatter(corr_score(:,1), corr_score(:,2))

[corr_coeff,corr_score,corr_latent,~,~,corr_mu] = pca(corr_fr);

figure; plot(corr_latent)




