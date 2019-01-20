%% Prompting whether to reload data 
choice_list = {'Load data from beginning and plot', 'Plot from already loaded data'}; 
fig_num_plt = 5; 
choices = questdlg(sprintf('Choices to plot Fig %d?', fig_num_plt), ...
    sprintf('Figure %d Plotting', fig_num_plt), ...
	choice_list{:}, choice_list{2});
% Handle response
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
    %% Loading raw data
    data_folder = '../data/extended_net_withGJ/';
    filename_prefix = 'extendedNet_withGJ_moredtvariations_0*.mat';
    
    thr_lat = 0.10; % of Src 
    cell_pref = {'Src', 'Int', 'Tgt'};
    [stat_dist, x_axis_actual, y_axis_actual] = returnStatDist(data_folder, filename_prefix, cell_pref, thr_lat);
end

%% Load layout
file_name  = 'Fig5_layout_official.svg';
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
    'VerticalAlignment', 'bottom', 'FontSize', 16, 'FontWeight', 'normal'};

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:});

create_text = @(tag_name, string_, extra_style) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, text_style{:}, extra_style{:});

create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);

gelec_cmap_factor = 0.8; 
%% Value vectors
t_init_std = x_axis_actual.vec;
g_elec_vec = y_axis_actual.vec;

name_x_ = '\sigma_{inp}'; 

%% Fig 5A
create_ann('ann_A', 'A');

%% Text for Fig 5B
add_text_style = {'fontweight', 'normal', 'fontsize', 16}; 
create_text('text_B@1', 'Int', add_text_style); 
create_text('text_B@2', 'Tgt', add_text_style); 

%% Fig 5B
f_limx = @(x) [-0.5,1]*4*x + [0,5];
f_limx_inset = @(x) [-2.6*x + 3, -2.6*x + 10]; 
f_limy = @(y) [-0.1, 2.5];
create_ann('ann_B', 'B');


std_tinit2plt = [1, 3, 5, 10];
cells2plt = {'Int', 'Tgt'};

nline_plt = length(g_elec_vec);
cmap = parula(nline_plt) * gelec_cmap_factor;

value_demo = struct(); 
value_demo.C = zeros(length(std_tinit2plt), nline_plt); 
value_demo.D = zeros(length(std_tinit2plt), nline_plt); 
value_demo.E = zeros(length(std_tinit2plt), nline_plt); 

for i = 1:length(std_tinit2plt)
    loc_i = find(abs(t_init_std - std_tinit2plt(i)) < eps);
    for j = 1:length(cells2plt)
        %% plotting distributions
        tag_ij = ['fig_B' num2str(i) '_' num2str(j)];
        ax_ij = create_fig(tag_ij);
        hold(ax_ij, 'on');
        
        loc_j = find(strcmp(cell_pref, cells2plt{j}));
        loc_src = find(strcmp(cell_pref, 'Src'));
        src_centers = [];
        src_ptsh = [];
        for k = 1:nline_plt
            stat_dist_ikj = stat_dist(loc_i,k,loc_j);
            plot(ax_ij, stat_dist_ikj.centers, stat_dist_ikj.ptsh, ...
                'color', cmap(k,:), 'linewidth', 1.5);
            src_ptsh = [src_ptsh; stat_dist(loc_i,k,loc_src).ptsh]; %#ok<AGROW>
            
            %% for figure 5CDE
            if strcmp(cells2plt{j}, 'Tgt') 
                value_demo.C(i, k) = stat_dist_ikj.lat2Src;
                value_demo.D(i, k) = stat_dist_ikj.pernspk;                 
                value_demo.E(i, k) = stat_dist_ikj.mean; 
            end
        end
        if j == 1
            plot(ax_ij, stat_dist(loc_i,1,loc_src).centers,...
                mean(src_ptsh,1), '--k', 'linewidth', 1.5);
        end
        
        xlim(f_limx(std_tinit2plt(i)));
        ylim(f_limy([]));
        set(ax_ij, 'fontsize', 10);
        %% conditional labelling
        if j == 1
            title(ax_ij, [name_x_ ' = ' num2str(std_tinit2plt(i)) ' ms'], ...
                'fontweight', 'normal', 'fontsize', 14);
        end
        if j == length(cells2plt)
            xlabel(ax_ij, 'Time (ms)', 'fontsize', 12);
        else
            set(ax_ij, 'xcolor', 'none');           
        end
        if i == 1
            ylabel(ax_ij, 'Density (norm)', 'fontsize', 12);
            set(ax_ij, 'ytick', [0, 1, 2], 'yticklabel', [0, 1, 2]);
        else
            set(ax_ij, 'ytick', ''); set(ax_ij, 'ycolor', 'none'); 
        end 
        set(ax_ij, 'TickLength', [0.025,0.025]); 
        
        
        
        %% inset to zoom in Tgt latency
        if strcmp(cells2plt{j}, 'Tgt')
            ax_inset = create_fig(['inset_B' num2str(i) '_' num2str(j)]) ;
            hold(ax_inset, 'on');
            for k = 1:nline_plt
                stat_dist_ikj = stat_dist(loc_i,k,loc_j);
                plot(ax_inset, stat_dist_ikj.centers, stat_dist_ikj.ptsh, ...
                    'color', cmap(k,:), 'linewidth', 1);
            end
            gray_clr = [1,1,1]*0.6; 
            plot(ax_inset, [-1, 10]*50, [1,1]*thr_lat, '--', 'color', gray_clr);
            set(ax_inset, 'box', 'on', 'xtick', '','ytick','');
            xlim(ax_inset, f_limx_inset(std_tinit2plt(i)));
            ylim(ax_inset, [-0.01, 0.15]);
        end
    end
end

colormap(ax_ij, cmap);
clrbar = colorbar(ax_ij);
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'box', 'off', 'linewidth', 0.001, 'position', pos_clrbar, 'fontsize', 10);
caxis(ax_ij, [min(y_axis_actual.vec), max(y_axis_actual.vec)]);
title(clrbar, y_axis_actual.name, 'fontsize', 14);
%% Fig 5CDE
lim_x = [min(g_elec_vec), max(g_elec_vec)] + [-1,1]*0.2; 
lim_y = @(ly) ly + 0.15 * abs(diff(ly)) * [-1,1];
x_label_name = y_axis_actual.name; 
num_lines = length(std_tinit2plt); 
cmap = [0,0,0;
        190,190,190; 
        95,200,211;
        200,0,0]/255; 
plot_props = struct(); 
plot_props.C = struct('title', 'Tgt latency', 'ylabel', 'ms'); 
plot_props.D = struct('title', '% Tgt spikes', 'ylabel', '%'); 
plot_props.E = struct('title', 'Mean Tgt spike time', 'ylabel', 'ms'); 
 
for tag_ = fieldnames(plot_props)'
    tag = tag_{:}; 
    create_ann(['ann_' tag], tag); 
    ax = create_fig(['fig_' tag]); 
    hold(ax, 'on'); 
    set(ax,'fontsize', 10, 'ticklength', [0.02,0.02]);
    value_ax = value_demo.(tag); 
    for i = 1:num_lines
        clr_i = cmap(i,:); 
        plot(ax, g_elec_vec, value_ax(i,:), ...
            '-', 'Color', clr_i, 'LineWidth', 2); 
    end
    if ~strcmp(tag, 'D')
        xlabel(ax, x_label_name, 'fontsize', 12);
    else
        set(ax, 'xtick', '');
    end
    xlim(ax,lim_x);
    ylabel(ax, plot_props.(tag).ylabel, 'fontsize', 12); 
    current_ylim = ylim(ax); 
    ylim(ax, lim_y(current_ylim)); 
    title(ax, plot_props.(tag).title, 'fontsize', 14, 'fontweight', 'normal');

end
    
colormap(ax, cmap);
clrbar = colorbar(ax);
pos_clrbar = layout_map('colorbar_CDE').normz_pos;
set(clrbar, 'box', 'off', 'linewidth', 0.001, 'position', pos_clrbar, ...
    'fontsize', 10, 'ticks', linspace(0,1,num_lines), 'ticklabels', std_tinit2plt); 
title(clrbar, x_axis_actual.name, 'fontsize', 12);