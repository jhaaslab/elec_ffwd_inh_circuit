%% Prompting whether to reload data
choice_list = {'Analyze data from beginning and plot', 'Plot from saved data'};
fig_num_plt = 8;
choices = questdlg(sprintf('Choices to plot Fig %d?', fig_num_plt), ...
    sprintf('Figure %d Plotting', fig_num_plt), ...
    choice_list{:}, choice_list{2});

switch choices
    case choice_list{1}
        choice = 'ananalyze_first';
    case choice_list{2}
        choice = 'loadplot_now';
end

%% Clear data + add necessary paths
switch choice
    case 'ananalyze_first'
        clc; clearvars -except choice; close all;
        addpath(genpath('figures'), genpath('functions'));
    case 'loadplot_now'
        clc; close all;
        addpath(genpath('figures'), genpath('functions'));
end
%% Loading depending on choice2plt
folder_prefix_pairs = {
    {'data/extended_net_withGJ/','extendedNet_withGJ_moredtvariations_*.mat'}, ...
    {'data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=1nS_*.mat'}, ...
    {'data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=3nS_*.mat'}, ...
    {'data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=5nS_*.mat'}};
saved_folder = 'data/info_anly';
num_cases = length(folder_prefix_pairs); 
stat_dist = cell(1,num_cases);
cell_pref = {'Src', 'Int', 'Tgt'};

common_file_pref = 'extra_inf_'; 
if strcmp(choice, 'ananalyze_first')
    %% Analyzing the data then saving them
    saved_folder = 'data/info_anly';
    bw_opt = 'fixed';
    smw_opt = 1; % meaning no smoothing, this works better 
    tic
    parpool(4)
    parfor i = 1 : num_cases
        data_folder = folder_prefix_pairs{i}{1};
        file_prefix = folder_prefix_pairs{i}{2};
        saved_name = fullfile(saved_folder, [common_file_pref file_prefix(1:(regexp(file_prefix, '(_\*.mat)')-1))]);
        dat_i = return_info_measures(data_folder, file_prefix, ...
            saved_name, cell_pref, bw_opt, smw_opt);
        stat_dist{i} = dat_i;
    end
    toc
else
    %% Loading the data
    for i = 1 : num_cases
        data_folder = folder_prefix_pairs{i}{1};
        file_prefix = folder_prefix_pairs{i}{2};
        saved_name = fullfile(saved_folder, [common_file_pref file_prefix(1:(regexp(file_prefix, '(_\*.mat)')-1))]);
        dat_i = load([saved_name '.mat']);
        stat_dist{i} = dat_i.info_res;
    end
end

%% Load layout
file_name  = 'Fig8_layout_v2.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor, ...
    'PaperUnits', 'centimeters','PaperPosition', [0, 0, width, height]*conv_factor); 

%% Annotation and axes styles
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 18, 'FontWeight', 'bold'};
text_style = {'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight', 'normal'};

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:});

create_text = @(tag_name, string_, extra_style) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, text_style{:}, extra_style{:});

create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);

%% Fig 8A&B
create_ann('ann_A', 'A');
create_ann('ann_B', 'B');

field2plt = {'mi_src_out', 'mi_over_hout'};  

right_lbl_names = struct(...
    't_init_std', struct('name', '\sigma_{inp}', 'unit', 'ms'), ...
    'g_elec_int_to_int', struct('name', '\SigmaG_{elec}', 'unit', 'nS'), ...
    'mi_src_out', struct('name', 'Mutual information', 'unit', 'bits', ...
        'limit', struct('raw', [1.5, 11], 'gain', [-1.25, 0.7])), ...
    'mi_over_hsrc', struct('name', 'I(X,Y)/I(X)', 'unit', 'bits'), ...
    'mi_over_hout', struct('name', 'Transmission efficiency', 'unit', 'bits', ...
         'limit', struct('raw', [0.5, 1.1]*100, 'gain', [-0.45, -0.03])), ...
    'entropy_lat_out2src',  struct('name', 'Entropy of latency', 'unit', 'bits', ...
         'limit', struct('raw', [-0.1, 6], 'gain', [-0.45, -0.03])) ...
     ); 

idx_x = 3; 
idx_difflines = 2; 
x_values = stat_dist{1}.value_var{idx_x}; 
diff_line_vals = stat_dist{1}.value_var{idx_difflines}; 

x_props = right_lbl_names.(stat_dist{1}.var_names{idx_x}); 
line_props = right_lbl_names.(stat_dist{1}.var_names{idx_difflines}); 

cmap_factor = 3;
cmap = jet(length(diff_line_vals)*cmap_factor);

loc_pref2plt = strcmp(cell_pref, 'Tgt') ;

for i = 1:length(field2plt)
    field_i = field2plt{i};
    name_i = right_lbl_names.(field_i).name;
    unit_i = right_lbl_names.(field_i).unit;
    lim_i = right_lbl_names.(field_i).limit.raw; 
    for k = 1:num_cases
        ax_plt = create_fig(['fig_A' num2str(i) '_' num2str(k)]);
        hold(ax_plt, 'on');
        loc2plt = cell(1,3);
        loc2plt{end} = loc_pref2plt;
        inf_k = stat_dist{k}.info_measures; 
        for j = 1:length(diff_line_vals)
            loc2plt{idx_difflines-1} = j;
            loc2plt{idx_x-1} = ':';
            dat2splt = [inf_k(loc2plt{:}).(field_i)];
            if contains(field_i, 'over')
                dat2splt = 100*dat2splt;
            end
            plot(ax_plt, x_values, dat2splt, '.-', 'color', cmap(j*cmap_factor,:), 'linewidth', 1.5, 'MarkerSize', 20);            

        end
        xlim(ax_plt, [min(x_values), max(x_values)] + 0.2*[-1,1]); 
        ylim(ax_plt, lim_i); 
        set(ax_plt, 'fontsize', 10); 
        if i == 1
            ylabel(ax_plt, unit_i, 'fontsize', 12); 
            set(ax_plt, 'ytick', [0, 3, 6, 9]);
        else
            ylabel(ax_plt, '%', 'fontsize', 12);
        end
        if k == num_cases
            xlabel(ax_plt, [x_props.name ' (' x_props.unit ')'], 'fontsize', 12);
        end
        if k == 1
            title(ax_plt, name_i, 'fontsize', 15, 'fontweight', 'normal');
        end
        if k ~= num_cases 
            set(ax_plt, 'xtick', '');
        end
        set(ax_plt, 'ticklength', [1,1]*0.015);
    end
    
end

colormap(cmap);
clrbar = colorbar;
clrbar_show = diff_line_vals([1,ceil(end/2),end]); 
pos_clrbar = layout_map('colorbar_A').normz_pos;
set(clrbar, 'box', 'off', 'position', pos_clrbar, ...
    'ticks', linspace(0,1,length(clrbar_show)), 'ticklabels', clrbar_show); 
title(clrbar, [line_props.name ' (' line_props.unit ')'], 'fontsize', 12);

