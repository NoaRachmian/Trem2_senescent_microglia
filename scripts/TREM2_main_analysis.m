%% Identification of senescent, TREM2-expressing microglia in aging and Alzheimer’s disease model mouse brain
% This script contains most of the analyses required to reproduce the 
% results presented in Rachmian et al., 2024.
%
% Each analysis corresponds to a single panel in the paper and can be run
% separately from the rest of the analyses.
%
% Before running the sections of the main anlyses, 
% it is required to run the setup section which is responsible to 
% load the processed data and set the relevent directories pathways.

%% Setup section - Setting scripts and data directories pathways
clear;clc;
warning('off');

% define data and scripts paths
data_path = 'C:\Users\noarac.WISMAIN\Desktop\trem2_nature_neuro_scripts\data\';
scripts_path = 'C:\Users\noarac.WISMAIN\Desktop\trem2_nature_neuro_scripts\scripts\';
addpath(scripts_path) % add the scripts folder path into MATLAB paths to run scripts from a different directory

% load color map
load([scripts_path,'magma_colormap.mat'])

% load xFAD processed cytof data
xFAD_data_path = [data_path,'5Xfad.xlsx'];
[~,~,xFAD_data_raw] = xlsread(xFAD_data_path);
[xFAD_protiens_levels_mat,xFAD_mouse_genotype_ind,xFAD_unique_cell_types,xFAD_unique_protiens] = sort_cytof_data(xFAD_data_raw,'basic');

% load and sorting DMhTAU processed cytof data loading 
DMhTAU_data_path = [data_path,'DMhTAU.xlsx'];
[~,~,DMhTAU_data_raw] = xlsread(DMhTAU_data_path);
[DMhTAU_protiens_levels_mat,DMhTAU_mouse_genotype_ind,DMhTAU_unique_cell_types,DMhTAU_unique_protiens] = sort_cytof_data(DMhTAU_data_raw,'basic');

% load and sorting old mice processed cytof data loading 
old_mice_data_path = [data_path,'Old mice.xlsx'];
[~,~,old_mice_data_raw] = xlsread(old_mice_data_path);
[old_mice_protiens_levels_mat,old_mice_mouse_genotype_ind,old_mice_unique_cell_types,old_mice_unique_protiens] = sort_cytof_data(old_mice_data_raw,'basic');

% load and sorting TREM2 processed cytof data loading 
Trem2_data_path = [data_path,'Trem2.xlsx'];
[~,~,Trem2_data_raw] = xlsread(Trem2_data_path);
[Trem2_protiens_levels_mat,Trem2_mouse_genotype_ind,Trem2_unique_cell_types,Trem2_unique_protiens] = sort_cytof_data(Trem2_data_raw,'trem2');

% load and sorting TREM2 processed cytof data loading 
ABT_data_path = [data_path,'ABT.xlsx'];
[~,~,ABT_data_raw] = xlsread(ABT_data_path);
[ABT_protiens_levels_mat,ABT_mouse_genotype_ind,ABT_unique_cell_types,ABT_unique_protiens] = sort_cytof_data(ABT_data_raw,'trem2');

% load and sorting ApoE processed cytof data loading 
ApoE_data_path = [data_path,'ApoE.xlsx'];
[~,~,ApoE_data_raw] = xlsread(ApoE_data_path);
[ApoE_protiens_levels_mat,~,ApoE_unique_cell_types,ApoE_unique_protiens] = sort_cytof_data(ApoE_data_raw,'apoe');

%% Figure 1A - average cytof expression matrix across old mice experiment

mean_old_mice_protiens_levels_mat = mean(old_mice_protiens_levels_mat,3,'omitnan');
figure('units','normalized','position',[0.3 0.3 0.5 0.35])
imagesc(mean_old_mice_protiens_levels_mat,[0 6])
set(gca,'xtick',1:32','xticklabels',xFAD_unique_protiens,'ytick',1:9,'yticklabels',old_mice_unique_cell_types,...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
colormap(turbo)
colorbar

%% Figure 1D - average cytof expression matrix across xFAD mice experiment

mean_xFAD_protiens_levels_mat = mean(xFAD_protiens_levels_mat,3,'omitnan');
figure('units','normalized','position',[0.3 0.3 0.5 0.35])
imagesc(mean_xFAD_protiens_levels_mat)
set(gca,'xtick',1:32','xticklabels',xFAD_unique_protiens,'ytick',1:9,'yticklabels',xFAD_unique_cell_types,...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
colormap(turbo)
colorbar

%% Figure 1G - comparing average expression levels for sensescent microglia cells across different experiments

mean_xFAD_protiens_levels_mat = mean(xFAD_protiens_levels_mat,3,'omitnan');
mean_DMhTAU_protiens_levels_mat = mean(DMhTAU_protiens_levels_mat,3,'omitnan');
mean_old_mice_protiens_levels_mat = mean(old_mice_protiens_levels_mat,3,'omitnan');

mean_senescent_microglia_all_groups = [mean_xFAD_protiens_levels_mat(4,1:31);...
                                      mean_DMhTAU_protiens_levels_mat(8,:);...
                                      mean_old_mice_protiens_levels_mat(3,1:31)];

figure('units','normalized','position',[0.3 0.3 0.4 0.25])
imagesc(mean_senescent_microglia_all_groups,'AlphaData',~isnan(mean_senescent_microglia_all_groups))
hold on
plot(xlim,[1.5 1.5],'color','k','linewidth',1.5)
plot(xlim,[2.5 2.5],'color','k','linewidth',1.5)
set(gca,'xtick',1:32','xticklabels',xFAD_unique_protiens,'ytick',1:3,'yticklabels',{'5xFAD','DMhTAU','Old mice'},...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
colormap(turbo)
colorbar


%% Figure 1H - Similarity in expression pattern of different cell types across experiments

xlabels = {'BAM','CD11b-','CD11b- SnCs','Microglia SnCs','Resting microglia'};
mean_xFAD_protiens_levels_mat = mean(xFAD_protiens_levels_mat,3,'omitnan')';
mean_DMhTAU_protiens_levels_mat = mean(DMhTAU_protiens_levels_mat,3,'omitnan')';
mean_old_mice_protiens_levels_mat = mean(old_mice_protiens_levels_mat,3,'omitnan')';

xFAD_relevent_cell_types = mean_xFAD_protiens_levels_mat(1:31,[5,8,9,4]);
DMhTAU_relevent_cell_types = mean_DMhTAU_protiens_levels_mat(1:31,[1,2,3,8]);
old_mice_relevent_cell_types = mean_old_mice_protiens_levels_mat(1:31,[4,7,8,3]);

% similarity between cell types across groups
relevent_cell_types_all_groups = {xFAD_relevent_cell_types',DMhTAU_relevent_cell_types',old_mice_relevent_cell_types'};
xlabels = {'CD11b-','CD11b- SnCs','Microglia SnCs','BAM'};
group_names = {'5xFAD','DMhTAU','Old mice'};

figure
similarity_across_all_groups = [];
ind = 1;
for group1 = 1:3
    current_group1 = relevent_cell_types_all_groups{group1};
    for group2 = 1:3
        current_group2 = relevent_cell_types_all_groups{group2};
        if group1 < group2 
        similarity_across_groups = corr(current_group1',current_group2','rows','pairwise');
        similarity_across_all_groups(:,:,ind) = similarity_across_groups;
        subplot(1,3,ind)
        imagesc(similarity_across_groups,[-0.2 1])
        set(gca,'xtick',1:5,'ytick',1:5)
        ylabel(group_names{group1})
        xlabel(group_names{group2})
        axis square
        ind = ind + 1;
        end
    end
end
colormap(magma_colormap)


%% Figure 2A - Expression levels of Trem2 in DAM and SnCs microglia

xFAD_TREM2_DAM = squeeze(xFAD_protiens_levels_mat(3,5,:));
xFAD_TREM2_SnCs = squeeze(xFAD_protiens_levels_mat(4,5,:));

xFAD_TREM2_both = [xFAD_TREM2_DAM,xFAD_TREM2_SnCs];

rng(7)
figure
turbo_colormap = colormap(turbo);
hold on
figure_boxplot(xFAD_TREM2_both)
for mouse = 1:size(xFAD_TREM2_both,1)
    x_jitter = randn(1)*0.05;
    plt(mouse) = scatter([1,2]+x_jitter,xFAD_TREM2_both(mouse,:),60,turbo_colormap(1+30*(mouse-1),:),'filled');
end
lgd = legend(plt,num2str([1:8]'),'Location','southeast');
title(lgd,'Mouse ID')
legend('boxoff')
ylim([0 2.75])
xlim([0.5 2.5])
axis square
set(gca,'xtick',1:2,'xticklabels',{'Disease-associated microglia','Senescent microglia'})
ylabel('Median expression (a.u)')

[~,p,~,stats] = ttest(xFAD_TREM2_DAM,xFAD_TREM2_SnCs)



%% Figure 2B - the relationship between the expression of TREM2 and different SnCs markers

[xFAD_trem2_sncs_mat] = trem2_sncs(xFAD_protiens_levels_mat);
[DMhTAU_trem2_sncs_mat] = trem2_sncs(DMhTAU_protiens_levels_mat);
[old_mice_trem2_sncs_mat] = trem2_sncs(old_mice_protiens_levels_mat);


exp_labels = [ones(1,size(xFAD_trem2_sncs_mat,1)),...
ones(1,size(DMhTAU_trem2_sncs_mat,1))*2,...
ones(1,size(old_mice_trem2_sncs_mat,1))*3];

comp_labels = {'p16','p19','p21'};

TREM2_across_exp = [xFAD_trem2_sncs_mat(:,1);DMhTAU_trem2_sncs_mat(:,1);old_mice_trem2_sncs_mat(:,1)];

figure('units','normalized','position',[0.3 0.3 0.6 0.25])
for comp = 1:3
    comp_gene_across_exp = [xFAD_trem2_sncs_mat(:,comp+1);DMhTAU_trem2_sncs_mat(:,comp+1);old_mice_trem2_sncs_mat(:,comp+1)];
    
    subplot(1,3,comp)
    
    hold on
    scatter(TREM2_across_exp,comp_gene_across_exp,[],exp_labels,'filled')
    h=lsline;
    h.LineWidth = 1.5;
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    [r,p] = corr(TREM2_across_exp,comp_gene_across_exp);
    text(-2.25,2.25,['r = ',num2str(round(r,3))])
    text(-2.25,2,['p = ',num2str(round(p,5))])
    ylabel('TREM2 expression (zscored)')
    xlabel([comp_labels{comp},' expression (zscored)'])
    axis square
end


%% Figure 2C - Expression profile for Trem2-/- and Trem2+/+ mice - only 5xFAD mice

Trem2_WT_ind = Trem2_mouse_genotype_ind == 3;
Trem2_KO_ind = Trem2_mouse_genotype_ind == 4;

Trem2_WT_mice_protiens_levels_mat = Trem2_protiens_levels_mat(:,:,Trem2_WT_ind);
Trem2_KO_mice_protiens_levels_mat = Trem2_protiens_levels_mat(:,:,Trem2_KO_ind);

mean_Trem2_WT_mice_protiens_levels_mat = mean(Trem2_WT_mice_protiens_levels_mat([1:4,6:7],:,:),3,'omitnan');
mean_Trem2_KO_mice_protiens_levels_mat = mean(Trem2_KO_mice_protiens_levels_mat([1:4,6:7],:,:),3,'omitnan');

figure('units','normalized','position',[0.05 0.3 0.8 0.3])
subplot(1,2,1)
imagesc(mean_Trem2_WT_mice_protiens_levels_mat,[0 6])
set(gca,'xtick',1:32','xticklabels',Trem2_unique_protiens,'ytick',1:6,'yticklabels',Trem2_unique_cell_types([1:4,6:7]),...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
title('5xFAD Trem2^+^/^+')

subplot(1,2,2)
imagesc(mean_Trem2_KO_mice_protiens_levels_mat,[0 6])
set(gca,'xtick',1:32','xticklabels',Trem2_unique_protiens,'ytick',[],...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
colormap(turbo)
colorbar
title('5xFAD Trem2^-^/^-')

%% Figure 5B - average expression levels across ABT mice

mean_ABT_protiens_levels_mat = mean(ABT_protiens_levels_mat,3,'omitnan');

figure('units','normalized','position',[0.3 0.3 0.5 0.35])
imagesc(mean_ABT_protiens_levels_mat,[0 6])
set(gca,'xtick',1:32','xticklabels',ABT_unique_protiens,'ytick',1:7,'yticklabels',ABT_unique_cell_types,...
    'TickLabelInterpreter', 'none','box','off')
xtickangle(90)
colormap(turbo)
colorbar


%% Extended data fig. 1 - relationship between ApoE, TREM2 and SnCs markers


[xFAD_trem2_sncs_mat] = trem2_sncs(xFAD_protiens_levels_mat);
[DMhTAU_trem2_sncs_mat] = trem2_sncs(DMhTAU_protiens_levels_mat);
[old_mice_trem2_sncs_mat] = trem2_sncs(old_mice_protiens_levels_mat);


exp_labels = [ones(1,size(xFAD_trem2_sncs_mat,1)),...
ones(1,size(DMhTAU_trem2_sncs_mat,1))*2,...
ones(1,size(old_mice_trem2_sncs_mat,1))*3];

comp_labels = {'TREM2','p16','p19','p21','p53','gH2AX'};

ApoE_across_exp = [xFAD_trem2_sncs_mat(:,7);DMhTAU_trem2_sncs_mat(:,7);old_mice_trem2_sncs_mat(:,7)];

figure('units','normalized','position',[0.3 0.3 0.8 0.45])
for comp = 1:4
    comp_gene_across_exp = [xFAD_trem2_sncs_mat(:,comp);DMhTAU_trem2_sncs_mat(:,comp);old_mice_trem2_sncs_mat(:,comp)];
  
    subplot(1,4,comp)
    hold on
    scatter(ApoE_across_exp,comp_gene_across_exp,[],exp_labels,'filled')
    h=lsline;
    h.LineWidth = 1.5;
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    [r,p] = corr(ApoE_across_exp,comp_gene_across_exp);
    text(-2.25,2.25,['r = ',num2str(round(r,3))])
    text(-2.25,2,['p = ',num2str(round(p,5))])
    ylabel('ApoE expression (zscored)')
    xlabel([comp_labels{comp},' expression (zscored)'])
    axis square
end

%% Extended data fig. 2 - Volcano plot between SnCs microglia and DAM

DAM_microglia = squeeze(xFAD_protiens_levels_mat(3,:,:));
SnCs_microglia = squeeze(xFAD_protiens_levels_mat(4,:,:));

exp_fold_change = [];
pval = [];
for prot = 1:size(DAM_microglia,1)
    current_prot_DAM = DAM_microglia(prot,:)';
    current_prot_SnCs = SnCs_microglia(prot,:)';
    
    avg_prot_DAM = mean(current_prot_DAM,'omitnan');
    avg_prot_SnCs = mean(current_prot_SnCs,'omitnan');
    
    exp_fold_change(prot) = avg_prot_SnCs./avg_prot_DAM;
    [~,pval(prot)] = ttest(current_prot_DAM,current_prot_SnCs);
end

corrected_alpha = 0.05./size(DAM_microglia,1);

sig_prot = exp_fold_change<2 & pval <= corrected_alpha;
sig_prot_high = exp_fold_change>=2 & pval <= corrected_alpha;
non_sig = ~(sig_prot|sig_prot_high);

all_sig = find(pval <= corrected_alpha);
chosen_prot = [all_sig(2:6),all_sig(end-4:end)];

figure('units','normalized','position',[0.3 0.3 0.5 0.5])
hold on
  xlim([-2 9])
scatter(log2(exp_fold_change(non_sig)),-log10(pval(non_sig)),40,[0.7 0.7 0.7],'filled')
scatter(log2(exp_fold_change(sig_prot_high)),-log10(pval(sig_prot_high)),40,[0.1 0.8 1],'filled')
scatter(log2(exp_fold_change(sig_prot)),-log10(pval(sig_prot)),40,[1 0.8 0.1],'filled')
plot([-1 -1],ylim,'--','color',[0.6 0.6 0.6])
plot([1 1],ylim,'--','color',[0.6 0.6 0.6])
plot(xlim,-log10([corrected_alpha corrected_alpha]),'--','color',[0.6 0.6 0.6])
ylabel('-log_1_0(p-value)')
xlabel('log_2(fold change)')

for prot = 1:length(chosen_prot)
    text(log2(exp_fold_change(chosen_prot(prot)))+0.025,-log10(pval(chosen_prot(prot)))+0.1,xFAD_unique_protiens{chosen_prot(prot)}(7:end),...
        'Interpreter', 'none')
end
