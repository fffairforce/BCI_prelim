function [output_data, PeakInfo, data_info, initPeak_flag_cat, X_F, X_F_cat, date_label] = DataReading(dates)
%reads neural and kinematic and lfads_output data  
%Input:     
% dates:            string to specify which datesto read in the data
%Output: 
% output_data:      lfads factor and rates extracted from LFADS model
% PeakInfo:         structure of neural data on speed peak on each input dates
% data_info:        structure of kinematic data on each input dates
% initPeak_flag_cat:initial trial flag on all input dates
% X_F:              lfads factors by day(tr F t)
% X_F_cat:          neural data
% date_label:
import PBT_analysis.*
addpath('D:\documents\2020 summer-LFAD\COT_LFAD\+PBT_analysis\')
for d=1:length(dates)
output_data=[];
lfads_output_dir = ['E:\Doc\COT_Data\multisession_lfads_output_0806\multisession_lfads_output_2017' dates{d}];
if ~exist(lfads_output_dir, 'dir')
    lfads_output_dir = ['D:\documents\Rouse Lab\GaTech colab\multisession_lfads_output_0806\multisession_lfads_output_' dates{d}];
    if ~exist(lfads_output_dir, 'dir')
        lfads_output_dir = ['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_LFADS_results\multisession_lfads_output_0806\multisession_lfads_output_' dates{d}];
    end
end
lfads_input_file = [lfads_output_dir '\lfads_2017' dates{d} '.h5'];%how to specify in multi-session
output_data = read_lfads_output(lfads_output_dir, lfads_input_file);
X_F{d}=permute(output_data.factors,[3 1 2]);%tr F t
lfads_input_file=[];
% X_F{d}=squeeze(X_F(:,:,16));
end
initPeak_flag_cat=[];
for d=1:length(dates)
PeakInfo_dir = ['E:\Doc\COT_Data\Data'];
data_info_dir = ['E:\Doc\COT_Data\PSTH'];
if ~exist(PeakInfo_dir, 'dir')
PeakInfo_dir = ['D:\documents\2020 summer-LFAD\COT_task\DataFiles\monk_p'];
end
if ~exist(PeakInfo_dir, 'dir')
PeakInfo_dir = ['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_LFADS_results\monk_p'];
end

if ~exist(data_info_dir, 'dir')
data_info_dir = ['D:\documents\2020 summer-LFAD\COT_task\DataFiles\monk_p\PSTH'];  
end
if ~exist(data_info_dir, 'dir')
data_info_dir = ['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_LFADS_results\PSTH'];  
end

PeakInfo{d} = load([PeakInfo_dir '\P_Spikes_2017' dates{d} '-data.mat'],'PeakInfo');
data_info{d} = load([data_info_dir '\P_Spikes_2017' dates{d} '-data_PSTH_prep_bin_20.mat']);
ind2include = squeeze(~isnan(data_info{d}.spikes_peakVel(1,1,:)));% remove nan trials 
initPeak_flag_cat = cat(1,initPeak_flag_cat,PeakInfo{d}.PeakInfo.initPeak_flag(ind2include)); 
end
% rng(33);
% n_corr = sum(initPeak_flag_cat==0);
% n_init_orig = sum(initPeak_flag_cat==1);
% init_trials = find(initPeak_flag_cat==1);
% corr_trials = find(initPeak_flag_cat==0);
% init_trials = init_trials(randperm(n_init_orig,n_corr));
% n_init = n_corr;
% all_trials = sort([init_trials;corr_trials]);
X_F_cat=[];
date_label=[];
for d=1:length(dates)
   X_F_tmp = X_F{d};
%    X_F_tmp = squeeze(X_F_tmp(:,:,16));
   X_F_cat = cat(1,X_F_cat,X_F_tmp);
   date_label = cat(1,date_label,d*ones(size(X_F_tmp,1),1));%maybe better way to do so?
end
