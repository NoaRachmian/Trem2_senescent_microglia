function [sorted_protiens_levels_mat,mouse_genotype_ind,unique_cell_types,unique_protiens] = sort_cytof_data(raw_xls_data,experiment_type)
    
    % extracting basic features from data
    data_raw_no_header = raw_xls_data(2:end,:);
    unique_mice_names = unique(data_raw_no_header(:,1));
    unique_cell_types = unique(data_raw_no_header(:,2),'stable');
    unique_protiens = unique(data_raw_no_header(:,3),'stable');

    sorted_protiens_levels_mat = nan(length(unique_cell_types),length(unique_protiens),length(unique_mice_names)); % 3D mat of cell typs x protiens x mice
    mouse_genotype_ind = [];
    for mouse_id = 1:length(unique_mice_names)
        current_mouse_name = unique_mice_names(mouse_id);
        current_mouse_rows = ismember(data_raw_no_header(:,1),current_mouse_name);
        current_mouse_data = data_raw_no_header(current_mouse_rows,:);

        if strcmp(experiment_type,'basic')
            
            if ~isempty(strfind(current_mouse_name{1},'WT'))
                mouse_genotype_ind(mouse_id) = 1; % is a WT
            elseif ~isempty(strfind(current_mouse_name{1},'AD'))
                mouse_genotype_ind(mouse_id) = 2; % is a AD
            end
            
        elseif strcmp(experiment_type,'trem2')
            
            if ~isempty(strfind(current_mouse_name{1},'WT_WT'))
                mouse_genotype_ind(mouse_id) = 1; % is a WT
            elseif ~isempty(strfind(current_mouse_name{1},'WT_KO'))
                mouse_genotype_ind(mouse_id) = 2; % is a WT KO
            elseif ~isempty(strfind(current_mouse_name{1},'AD_WT'))
                mouse_genotype_ind(mouse_id) = 3; % is a AD WT
            elseif ~isempty(strfind(current_mouse_name{1},'AD_KO'))
                mouse_genotype_ind(mouse_id) = 4; % is a AD KO
            end
            
        elseif strcmp(experiment_type,'ABT')
            
            if ~isempty(strfind(current_mouse_name{1},'WT'))
                mouse_genotype_ind(mouse_id) = 1; % is a WT
            elseif ~isempty(strfind(current_mouse_name{1},'DMSO'))
                mouse_genotype_ind(mouse_id) = 2; % is a DMSO
            elseif ~isempty(strfind(current_mouse_name{1},'ABT'))
                mouse_genotype_ind(mouse_id) = 3; % is a ABT
            end
            
        end
        
        
            for cell_type_id = 1:length(unique_cell_types)
            current_cell_type_name = unique_cell_types(cell_type_id);
            current_cell_type_rows = ismember(current_mouse_data(:,2),current_cell_type_name);
            current_cell_type_data = current_mouse_data(current_cell_type_rows,:);

            sorted_protiens_levels_mat(cell_type_id,:,mouse_id) = cell2mat(current_cell_type_data(:,end));
        end

    end

end