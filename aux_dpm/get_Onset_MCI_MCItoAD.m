function [ages_onset,gt_onset,metrics]=get_Onset_MCI_MCItoAD(DPS_test,posterior_bayes,ages_test,labels_test,btstrp)
% ============================================================
% Project:    Disease progression modeling from MCI to Dementia
% Repository: https://github.com/cplatero/MCItoDementia_DPM
% Author:     Carlos Platero
% Email:      carlos.platero@upm.es
% Institution:Universidad Politécnica de Madrid 
% ------------------------------------------------------------
% Filename:    get_Onset_MCI_MCItoAD.m
% Description: Convert time
% 
% Version:    1.0
% Date:       2025-05-09
% MATLAB Ver: R2024a 
% ============================================================
I=size(DPS_test,1);
ages_onset=nan(I,btstrp);
for n = 1 : btstrp
    dps_test_n = DPS_test(:, :, n);
    prob_MCI=posterior_bayes{2, n}(dps_test_n)>.5;
    for i=1: I
        if(sum(prob_MCI(i,:)>0))
            mask_values= ~isnan(ages_test(i,:));
            age_subj=ages_test(i,mask_values);
            label_subj=prob_MCI(i,mask_values);
            idx=find(label_subj);
            if(label_subj(end))
                if(idx(1)>1)
                    ages_onset(i,n)=(age_subj(idx(1))+age_subj(idx(1)-1))/2;
                else
                    ages_onset(i,n)=age_subj(idx(1));
                end
            end
        end %convert
    end %subjects
end %btstrp

gt_onset=get_GT_onset(labels_test,ages_test);

%% metrics
% sMCI + AD
age_bsl=ages_test(:,1);

metrics=table;
percentage_sMCI=nan(btstrp,1);
percentage_AD=nan(btstrp,1);
percentage_pMCI=nan(btstrp,1);
corr_MCI_age=nan(btstrp,1);
corr_MCI_reserve=nan(btstrp,1);

num_sMCI=sum(isnan(gt_onset));
num_AD=sum(gt_onset==age_bsl);

for n = 1 : btstrp
    percentage_sMCI(n)=sum(isnan(ages_onset(:,n)) & isnan(gt_onset))/num_sMCI*100;
    percentage_AD(n)=sum((ages_onset(:,n)==age_bsl) & (gt_onset==age_bsl))/num_AD*100;
end

% pMCI
num_pMCI=sum(gt_onset>age_bsl);

for n = 1 : btstrp
    percentage_pMCI(n)=sum((ages_onset(:,n)>age_bsl) & (gt_onset>age_bsl))/num_pMCI*100;
    cov=corrcoef(rmmissing([ages_onset(:,n),gt_onset]));
    try
        corr_MCI_age(n)=cov(1,2);
    catch
        corr_MCI_age(n)=0;
    end

    cov=corrcoef(rmmissing([ages_onset(:,n)-age_bsl,gt_onset-age_bsl]));
    try
        corr_MCI_reserve(n)=cov(1,2);
    catch
        corr_MCI_reserve(n)=0;
    end
    
end

metrics.percentage_sMCI=percentage_sMCI;
metrics.percentage_AD=percentage_AD;
metrics.percentage_pMCI=percentage_pMCI;
metrics.corr_MCI_age=corr_MCI_age;
metrics.corr_MCI_reserve=corr_MCI_reserve;

%% show
fprintf('Percentage sMCI: %.1f (%.2f)\n',mean(metrics.percentage_sMCI),...
    std(metrics.percentage_sMCI));
% fprintf('Percentage AD: %.1f (%.2f)\n',mean(metrics.percentage_AD),...
%     std(metrics.percentage_AD));
fprintf('Percentage pMCI: %.1f (%.2f)\n',mean(metrics.percentage_pMCI),...
    std(metrics.percentage_pMCI));

fprintf('Correlation age conversion: %.2f (%.2f)\n',mean(metrics.corr_MCI_age),...
    std(metrics.corr_MCI_age));
fprintf('Correlation cognition reserve: %.2f (%.2f)\n',mean(metrics.corr_MCI_reserve),...
    std(metrics.corr_MCI_reserve));


end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gt_onset=get_GT_onset(labels_test,ages_test)
I=size(labels_test,1);
gt_onset=nan(I,1);
% MCI=strcmp(labels_test,'Dementia');
MCI=strcmp(labels_test,'SI');

for i=1: I
    if(sum(MCI(i,:)>0))
        mask_values= ~isnan(ages_test(i,:));
        age_subj=ages_test(i,mask_values);
        label_subj=MCI(i,mask_values);
        idx=find(label_subj);
        if(label_subj(end))
            if(idx(1)>1)
                gt_onset(i)=(age_subj(idx(1))+age_subj(idx(1)-1))/2;
            else
                gt_onset(i)=age_subj(idx(1));
            end
        end
    end %convert
end %subjects

end