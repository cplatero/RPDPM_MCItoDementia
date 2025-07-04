function multi_RPDPM_MCItoDementia
% ============================================================
% Project:    Disease progression modeling from MCI to Dementia
% Repository: https://github.com/cplatero/MCItoDementia_DPM
% Author:     Carlos Platero
% Email:      carlos.platero@upm.es
% Institution:Universidad Politécnica de Madrid 
% ------------------------------------------------------------
% Filename:    multi_RPDPM_MCItoDementia.m
% Description: Script for training the RPDPM model 
%              on longitudinal data.
% 
% Version:    1.0
% Date:       2025-05-09
% MATLAB Ver: R2024a 
% ============================================================
clc;
close all;
addpath('./rpdpm');
addpath('./aux_dpm');

%% Input data file


% Markers:
    % 01 - IB (-)
    % 02 - IL (-) 
    % 03 - MEC (-)
    % 04 - mec30(-)
    % 05 - FV (-)
    % 06 - Reloj (-)
    % 07 - EURO_D (-)
    % 08 "mir342"
    % 09 "mir128"
    % 10 "mir132"
    % 11 "mir146"
    % 12 "mir155"
    % 13 "mir181"
    % 14 "mir206"
    % 15 "mir103"
    % 16 "mir7"
    % 17 "ptau"
    % 18 "NFL"
    % 19 "gfap"
    % 20 "tautot"
    % 21 "AB40"
    % 22 "AB42"
    % 23 "Abratio"

%vector_range = 2:5;   
vector_range = 3;
features_set = 1:6;
head_fich_name ="./output/MCItoDementia__";
%btstrp=10;
btstrp=2;
classes={'NO','SI'};

% Ranges of biomarkers
ranges = [...
    0 100;... IB
    0 4;... IL
    0 35;... MEC
    0 30;... mec30
    0 19;...FV
    0 9;...Reloj
    0 9;...EURO_D
    -1.4854    6.5275;..."mir342"
     0.6428    7.3095;..."mir128"
    -0.3805   10.5409;..."mir132"
    -2.2297    4.3806;..."mir146"
    -1.0223   12.0812;..."mir155"
    -10.000   10.9117;..."mir181"
    -5.0000   19.5920;..."mir206"
    -3.7111    4.3548;..."mir103"
    -0.9568    6.6876;..."mir7"
     5.6755   72.5706;..."ptau"
     4.7295   53.6308;..."NFL"
     25.4255  400.00;..."gfap"
     0.9209   10.1192;..."tautot"
     44.1972  400.000;..."AB40"
     1.0382   12.4139;..."AB42"
     0.0113    0.0646;...."Abratio"
    ];


load('./data/data_MCItoAD_RNA_model_multi_full','y','ages','labels','diagn',...
    'idx_train','idx_test','feature_names');


%% Pool marker sets    
features_set=1:length(features_set);
feature_names=feature_names(features_set);
max_dim=vector_range(end);
features_subset_matrix = nan(10000,max_dim);
next_row=1;
for i = 1 : length(vector_range)
    aux_subsets = nchoosek(features_set,vector_range(i));
    features_subset_aux=nan(size(aux_subsets,1),max_dim);
    features_subset_aux(:,1:vector_range(i))=aux_subsets;
    n_rows=size(features_subset_aux,1);
    features_subset_matrix(next_row:n_rows+next_row-1,:)=features_subset_aux;
    next_row=next_row+n_rows;
end
features_subset_matrix(next_row:end,:)=[];

%% Read models
alreadyBuiltDPM = false;
if alreadyBuiltDPM
    mask_DPM_built=false(size(features_subset_matrix,1),1);
    for i=1:size(features_subset_matrix,1)
        idx_markers=features_subset_matrix(i,:);
        idx_markers=idx_markers(~isnan(idx_markers));
        str_feat=strrep(int2str(idx_markers), ' ', '_');
        fich_name=head_fich_name+str_feat;
        if(exist(fich_name+'.mat','file'))
            mask_DPM_built(i)=true;
        end
    end
    features_subset_matrix=features_subset_matrix(mask_DPM_built,:);
end

%% Performance models
mean_MAE=nan(size(features_subset_matrix,1),length(features_set),2);
std_MAE= nan(size(features_subset_matrix,1),length(features_set),2);
mean_NMAE=nan(size(features_subset_matrix,1),2);
std_NMAE= nan(size(features_subset_matrix,1),2); 
mean_AUC=nan(size(features_subset_matrix,1),2);
std_AUC= nan(size(features_subset_matrix,1),2);

mean_percent_sMCI=nan(size(features_subset_matrix,1),2);
std_percent_sMCI=nan(size(features_subset_matrix,1),2);
mean_percent_AD=nan(size(features_subset_matrix,1),2);
std_percent_AD=nan(size(features_subset_matrix,1),2);
mean_percent_pMCI=nan(size(features_subset_matrix,1),2);
std_percent_pMCI=nan(size(features_subset_matrix,1),2);

mean_corr_MCI_age=nan(size(features_subset_matrix,1),2);
std_corr_MCI_age=nan(size(features_subset_matrix,1),2);
mean_corr_MCI_reserve=nan(size(features_subset_matrix,1),2);
std_corr_MCI_reserve=nan(size(features_subset_matrix,1),2);

str_feat=string.empty(size(features_subset_matrix,1), 0);

%% Building DPMs
% pc = parcluster('local');
% pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
% parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

for i = 1 : size(features_subset_matrix,1)
    idx_markers=features_subset_matrix(i,:);
    idx_markers=idx_markers(~isnan(idx_markers));
    model_feat=strjoin(feature_names(idx_markers), ', ');
    fprintf('Model %d of %s\n', i, model_feat);
    str_feat(i)=strrep(int2str(idx_markers), ' ', '_');
    fich_name=head_fich_name+str_feat(i);
    if(~exist(fich_name+'.mat','file'))
        [MAE,NMAE,AUC,percent_sMCI,percent_AD,percent_pMCI,corr_MCI_age,corr_MCI_reserve]=...
            build_rpdpm_MCItoAD(y(idx_markers,:,:),ages,labels,diagn,...
            classes,ranges(idx_markers,:),idx_train,idx_test,btstrp,...
            fich_name,feature_names(idx_markers));
       
        
    else
        load(fich_name+'.mat', 'auc','MAE','NMAE','m_test','m_train');
        AUC=auc;
        percent_sMCI=[m_test.percentage_sMCI,m_train.percentage_sMCI];
        percent_AD=[m_test.percentage_AD,m_train.percentage_AD];
        percent_pMCI=[m_test.percentage_pMCI,m_train.percentage_pMCI];
        corr_MCI_age=[m_test.corr_MCI_age,m_train.corr_MCI_age];
        corr_MCI_reserve=[m_test.corr_MCI_reserve,m_train.corr_MCI_reserve];
        
    end

        
    mean_NMAE(i,:)=mean(NMAE);
    std_NMAE(i,:)=std(NMAE);
    mean_MAE(i,idx_markers,:)=mean(MAE);
    std_MAE(i,idx_markers,:)=std(MAE);
    mean_AUC(i,:)=mean(AUC);
    std_AUC(i,:)=std(AUC);

    mean_percent_sMCI(i,:)=mean(percent_sMCI);
    std_percent_sMCI(i,:)=std(percent_sMCI);
    mean_percent_AD(i,:)=mean(percent_AD);
    std_percent_AD(i,:)=std(percent_AD);
    mean_percent_pMCI(i,:)=mean(percent_pMCI);
    std_percent_pMCI(i,:)=std(percent_pMCI);


    mean_corr_MCI_age(i,:)=mean(corr_MCI_age);
    std_corr_MCI_age(i,:)=std(corr_MCI_age);
    mean_corr_MCI_reserve(i,:)=mean(corr_MCI_reserve);
    std_corr_MCI_reserve(i,:)=std(corr_MCI_reserve);
        
    fprintf('\n');
end
    
    

%% Overall results
results=table;
results.str_feat=str_feat';
results.mean_percent_sMCI=mean_percent_sMCI;
results.mean_percent_AD=mean_percent_AD;
results.mean_percent_pMCI=mean_percent_pMCI;
results.mean_corr_MCI_age=mean_corr_MCI_age;
results.mean_corr_MCI_reserve=mean_corr_MCI_reserve;
results.mean_AUC=mean_AUC;

save('Metrics_multi','mean_MAE','std_MAE','mean_NMAE','std_NMAE','std_AUC',...
    'std_percent_sMCI','std_percent_AD','std_percent_pMCI',...
    'std_corr_MCI_age','std_corr_MCI_reserve','results',...
    'feature_names','vector_range');







end



