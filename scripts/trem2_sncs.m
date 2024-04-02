function [trem2_sncs_mat] = trem2_sncs(protiens_levels_mat)

TREM2 = squeeze(protiens_levels_mat(:,5,:));
p16 = squeeze(protiens_levels_mat(:,10,:));
p19 = squeeze(protiens_levels_mat(:,11,:));
gH2AX = squeeze(protiens_levels_mat(:,13,:));
p21 = squeeze(protiens_levels_mat(:,14,:));
p53 = squeeze(protiens_levels_mat(:,15,:));
ApoE = squeeze(protiens_levels_mat(:,21,:));


norm_avg_TREM2 = zscore(nanmean(TREM2))';
norm_avg_p16 = zscore(nanmean(p16))';
norm_avg_p19 = zscore(nanmean(p19))';
norm_avg_p21 = zscore(nanmean(p21))';
norm_avg_p53 = zscore(nanmean(p53))';
norm_avg_gH2AX = zscore(nanmean(gH2AX))';
norm_avg_ApoE = zscore(nanmean(ApoE))';

trem2_sncs_mat = [norm_avg_TREM2,norm_avg_p16,norm_avg_p19,norm_avg_p21,norm_avg_p53,norm_avg_gH2AX,norm_avg_ApoE];
end