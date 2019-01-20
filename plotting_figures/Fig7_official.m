%% Prompting whether to reload data
choice_list = {'Load data from beginning and plot', 'Plot from already loaded data'};
fig_num_plt = 7;
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
        '../data/extended_net_withGJ/','extendedNet_withGJ_moredtvariations_0*.mat';
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=1nS_*.mat';
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=3nS_*.mat'; 
        '../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=5nS_*.mat'};
    num_cases = size(folder_prefix_pairs, 1);
    stat_dist = cell(1,num_cases);
        
    thr_lat = 0.10; % of Src 
    cell_pref = {'Src', 'Int', 'Tgt'};
    for i = 1 : num_cases
        data_folder = folder_prefix_pairs{i,1};
        file_prefix = folder_prefix_pairs{i,2};
        [stat_dist_i, x_axis_actual, y_axis_actual] = returnStatDist(data_folder, file_prefix, cell_pref, thr_lat);
        stat_dist{i} = stat_dist_i;
    end
end
%% Load layout
file_name  = 'Fig7_layout_official.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor, ...
    'PaperUnits', 'centimeters','PaperPosition', [0, 0, width, height]*conv_factor, ...
    'PaperSize', [width, height]*conv_factor);

%% Annotation and axes styles
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 20, 'FontWeight', 'bold'};
text_style = {'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'normal'};

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:});

create_text = @(tag_name, string_, extra_style) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, text_style{:}, extra_style{:});

create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);

%% Common colormap 
color_options = {
    [68,0,170]/255;
    [102,0,255]/255;
    [153,85,255]/255;
    [204,170,255]/255;
    [1,1,1]*0.999;
    [255,204,170]/255;
    [255,153,85]/255;
    [255,102,0]/255;
    [170,68,0]/255;
    };

cmap = createmap(color_options, 500, 0.2, 0.95);
pos_white = find(cmap(:,1)>0.99&cmap(:,2)>0.99&cmap(:,3)>0.99);
cmap(pos_white,:) = 0.99*cmap(pos_white,:);
%% Plot Fig 7 
flip_xy_axis = true;
caxis_mode = 'manual'; 

title_textstyle = {'fontsize', 15,  'fontweight', 'bold'}; 
create_ann('ann_A', 'A');
create_text('text_A', 'Spike time distribution of Int', title_textstyle);

create_ann('ann_B', 'B');
create_text('text_B', 'Spike time distribution of Tgt', title_textstyle);

stat_fields2plt = {'mean', 'std', 'peak', 'auc', 'lat2Src'};
right_lbl_names = struct(...
    'mean', struct('name', 'Mean spike time', 'unit', 'gain (ms)', ...
        'caxis_lim', struct('Int',[-6,2], 'Tgt', [-2.5,1.5])), ...
    'std', struct('name', 'Standard deviation', 'unit', 'gain (norm)', ...
        'caxis_lim', struct('Int',[-1.2,0.2], 'Tgt', [0,3])), ...
    'peak', struct('name', 'Max density', 'unit', 'gain (norm)', ...
        'caxis_lim', struct('Int',[-0.5,1.5], 'Tgt', [-1,0])), ...
    'auc', struct('name', 'Total response', 'unit', 'gain (norm)',...
        'caxis_lim',struct('Int',[-0.8,0.2], 'Tgt', [-0.4,0])), ...
    'lat2Src', struct('name', 'Latency', 'unit', 'gain (ms)', ...
        'caxis_lim', struct('Int', [0,20], 'Tgt', [-1.5,0]))...
    );
if flip_xy_axis
    x_axis_ = y_axis_actual;
    y_axis_ = x_axis_actual;
else
    x_axis_ = x_axis_actual;
    y_axis_ = y_axis_actual;
end
cells2plt = {'Int', 'Tgt'};
ann_cell = {'A', 'B'}; 
for i = 1:length(stat_fields2plt)
    field_i = stat_fields2plt{i};
    name_i = right_lbl_names.(field_i).name;
    unit_i = right_lbl_names.(field_i).unit;
    
    for j = 1:length(cells2plt)
        loc_j = find(strcmp(cell_pref, cells2plt{j}));
        
        for m = 1:num_cases
            
            data2plt_ij = arrayfun(@(x) x.(field_i), stat_dist{m}(:,:,loc_j));
            data2plt_ij = reshape(data2plt_ij, [size(data2plt_ij,1), size(data2plt_ij,2)]);
            
            data2plt_ij = data2plt_ij - repmat(data2plt_ij(:,1),[1,size(data2plt_ij,2)]);
            %% plotting heat map
            ax_pixels = create_fig(['fig_' ann_cell{j} num2str(i) '_' num2str(m)]);
            hold(ax_pixels, 'on');
            x_vec = x_axis_.vec;
            y_vec = y_axis_.vec;
            if flip_xy_axis
                data2plt_ij = data2plt_ij';
            end
            image(ax_pixels, x_vec, y_vec, data2plt_ij', 'CDataMapping', 'scaled');
            set(ax_pixels, 'ydir', 'normal', 'fontsize', 10);
            xlim(ax_pixels, [min(x_vec), max(x_vec)] + 0.52*[-1,1]*(x_vec(2)-x_vec(1)));
            ylim(ax_pixels, [min(y_vec), max(y_vec)] + 0.52*[-1,1]*(y_vec(2)-y_vec(1)));
            pbaspect(ax_pixels, [1,1,1]);
            %% colormap
            colormap(ax_pixels, cmap);
            
            if strcmp(caxis_mode, 'auto')
                caxis(ax_pixels, [min(data2plt_ij(:)), max(data2plt_ij(:))]);
            else
                caxis_ = right_lbl_names.(field_i).caxis_lim.(cells2plt{j});
                caxis_ = max(abs(caxis_));
                caxis(ax_pixels, caxis_*[-1,1] );
            end
            clrbar = colorbar; 
            set(clrbar,'visible','off');
            if m == num_cases
                clrbar = colorbar(ax_pixels);
                ylabel(clrbar, unit_i, 'fontsize', 12);
                set(clrbar, 'box', 'off');                
                
            end
            %% conditional labeling
            if i == 1 && m == num_cases
                set(ax_pixels, 'xtick', x_axis_.show([1,end]), 'xticklabel', x_axis_.show([1,end]));
                set(ax_pixels, 'ytick', y_axis_.show([1,end]), 'yticklabel', y_axis_.show([1,end]));
                xlabel(ax_pixels, x_axis_.name, 'fontsize', 14);
                ylabel(ax_pixels, y_axis_.name, 'fontsize', 14);
            else 
                 set(ax_pixels, 'xtick', '', 'ytick', '');  
            end
            if m == 1
                title(ax_pixels, name_i, 'fontsize', 12, 'fontweight', 'normal');
            end
            set(ax_pixels, 'linewidth', 1.5, 'box', 'on');
        end
    end
end

