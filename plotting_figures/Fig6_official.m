%% Prompting whether to reload data
choice_list = {'Load data from beginning and plot', 'Plot from already loaded data'};
fig_num_plt = 6;
choices = questdlg(sprintf('Choices to plot Fig %d?', fig_num_plt), ...
    sprintf('Figure %d Plotting', fig_num_plt), ...
    choice_list{:}, choice_list{2});

switch choices
    case 'Load data from beginning and plot'
        choice = 'loadfirst';
    case 'Plot from already loaded data'
        choice = 'plotnow';
end

%% Clear data + add necessary paths
switch choice
    case 'loadfirst'
        clc; clearvars -except choice; close all;
        addpath(genpath('../figures'), genpath('../functions'));
    case 'plotnow'
        clc; close all;
end

%% Loading and normalizing data depending on choice2plt
if strcmp(choice, 'loadfirst')
    %% Loading data
    folder_prefix_pairs = {
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=1nS_*.mat';
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=3nS_*.mat'; 
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=5nS_*.mat';};
    num_cases = size(folder_prefix_pairs, 1);
    stat_dist = cell(1,num_cases);
    cell_pref = {'Src', 'Int', 'Tgt'};
    for i = 1 : num_cases
        data_folder = folder_prefix_pairs{i,1};
        file_prefix = folder_prefix_pairs{i,2};
        [stat_dist_i, x_axis_actual, y_axis_actual] = returnStatDist(data_folder, file_prefix, cell_pref);
        stat_dist{i} = stat_dist_i;
    end
end

%% Load layout
file_name  = 'Fig6_layout_v4.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor, ...
    'PaperUnits', 'centimeters','PaperPosition', [0, 0, width, height]*conv_factor, 'PaperSize', [width, height]*conv_factor);

%% Annotation and axes styles
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 20, 'FontWeight', 'bold'};
text_style = {'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 18, 'FontWeight', 'normal'};

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:});

create_text = @(tag_name, string_, extra_style) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, text_style{:}, extra_style{:});

create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);

%% Value vectors 
x_axis_actual.name = '\sigma_{inp}';
x_axis_actual.unit = 'ms'; 

y_axis_actual.name = '\SigmaG_{elec}';
y_axis_actual.unit = 'nS'; 

t_init_std = x_axis_actual.vec; 
g_elec_vec = y_axis_actual.vec; 
g_gaba_int_vec = [1,3,5]; 

%% Fig 6A
create_ann('ann_A', 'A'); 

%% Plot Fig 6B
cell2plt = 'Tgt'; 
create_ann('ann_B', 'B'); 
create_text('text_B', [cell2plt ' spike time distributions'], {}); 
switch upper(cell2plt)
    case 'TGT'
        f_limx = @(x) [-1,0.8]*3.8*x + [5,10];
        f_limy = @(y) [-0.1, 1.75];
    case 'INT'
        f_limx = @(x) [-0.5,1]*4*x + [0,5];
        f_limy = @(y) [-0.1, 2.5];
end
std_tinit2plt = [1, 3, 5, 10];
loc_cell = find(strcmp(cell_pref, cell2plt));

nline_plt = length(g_elec_vec); 
cmap = parula(nline_plt); 

for i = 1:length(std_tinit2plt) 
    loc_i = find(abs(t_init_std - std_tinit2plt(i)) <eps);
    for j = 1:length(g_gaba_int_vec) 
        
        %% plotting the normalized distributions
        tag_ij = ['fig_B' num2str(i) '_' num2str(j)]; 
        ax_ij = create_fig(tag_ij); 
        hold(ax_ij, 'on'); 
        for k = 1:nline_plt
            stat_ijk = stat_dist{j}(loc_i,k,loc_cell);
            plot(ax_ij, stat_ijk.centers, stat_ijk.ptsh, ...
                'color', cmap(k,:), 'linewidth', 1);             
        end
        xlim(ax_ij, f_limx(std_tinit2plt(i))); 
        ylim(ax_ij, f_limy([]));
        set(ax_ij, 'fontsize', 10);
        
        %% conditional labeling
        if j == 1
            title(ax_ij, [x_axis_actual.name ' = ' num2str(std_tinit2plt(i)) ' ' x_axis_actual.unit], ...
                'fontweight', 'normal', 'fontsize', 14);
        end
        if j == length(g_gaba_int_vec)
            xlabel(ax_ij, 'Time (ms)', 'fontsize', 12);
        else
            set(ax_ij, 'xcolor', 'none');
        end
        if i == 1
            ylabel(ax_ij, 'Density (norm)', 'fontsize', 12);
            set(ax_ij, 'ytick', [0,1,2], 'yticklabel', [0,1,2]);
        else
            set(ax_ij, 'ycolor', 'none');
        end
        set(ax_ij, 'ticklength', [1,1]*0.025);
        
    end
end
colormap(ax_ij, cmap);
clrbar = colorbar(ax_ij);
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'box', 'off', 'linewidth', 0.001, 'position', pos_clrbar, 'fontsize', 10);
caxis(ax_ij, [min(y_axis_actual.vec), max(y_axis_actual.vec)]);
title(clrbar, [y_axis_actual.name ' (' y_axis_actual.unit ')'],...
    'fontsize', 14);
