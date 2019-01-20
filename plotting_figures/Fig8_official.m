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
        addpath(genpath('../figures'), genpath('../functions'));
    case 'loadplot_now'
        clc; close all;
        addpath(genpath('figures'), genpath('functions'));
end
%% Loading depending on choice2plt
folder_prefix_pairs = {
    {'../data/extended_net_withGJ/','extendedNet_withGJ_moredtvariations_*.mat'}, ...
    {'../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=1nS_*.mat'}, ...
    {'../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=3nS_*.mat'}, ...
    {'../data/extended_net_withGJandGABA/', 'extendedNet_withGJandGABA=5nS_*.mat'}};
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
file_name  = 'Fig8_layout_official.svg';
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

%% Some general labels 
right_lbl_names = struct(...
    't_init_std', struct('name', '\sigma_{inp}', 'unit', 'ms'), ...
    'g_elec_int_to_int', struct('name', '\SigmaG_{elec}', 'unit', 'nS'), ...
    'mi_src_out', struct('name', 'Mutual information', 'unit', 'bits', ...
        'limit', [3.7, 9.5]), ...
    'mi_over_hout', struct('name', 'Transmission efficiency', 'unit', '%', ...
         'limit', [50, 100])); 
     
%% Load demo for Fig8A 
demo_filepref_idx = 1; % without GABA
data_folder = folder_prefix_pairs{demo_filepref_idx}{1};
filename_prefix = folder_prefix_pairs{demo_filepref_idx}{2};
[raw_co_spike, num_unit_trial, sigma_inp, g_elec] = return_co_spike(data_folder, filename_prefix); 

%% Fig 8A 
create_ann('ann_A', 'A');
for i = 1:4
    create_text(['ann_A' num2str(i)], num2str(i), {'fontweight', 'bold'});
end

cellpop_lbl = struct();
cellpop_lbl.Src = struct('order', 1, 'color', [1,1,1]*0.8);
cellpop_lbl.Int = struct('order', 2, 'color', [0.1,0.75,0.2]);
cellpop_lbl.Tgt = struct('order', 3, 'color', [0.0,0.1,0.75]);
cell_keys = fieldnames(cellpop_lbl);

binw = 0.01;
edges = 0:binw:20;
centers = (edges(1:end-1)+edges(2:end))/2;

smooth_win = hanning(20);
smooth_win = smooth_win / sum(smooth_win);

raster_unit_trial = 1:num_unit_trial;
common_scatter_style = {'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none'};
marker_size = 20;

loc_cell_demo = [2,3]; % only int and tgt 
loc_sigma_demo = 5; % 5ms 

loc_gelec_demo = [1,3,5,7,10]; % gelec = 0, 1, 2, 3, 4.5 nS 
loc_gelec_raster = [3,10]; % gelec = 1, 4.5 nS 
n_gelec_demo = length(loc_gelec_demo); 

demo_raster_1 = create_fig('fig_A1'); hold(demo_raster_1, 'on'); 
demo_raster_2 = create_fig('fig_A2'); hold(demo_raster_2, 'on'); 
demo_raster = [demo_raster_1, demo_raster_2]; 
demo_distr  = create_fig('fig_A3');  hold(demo_distr, 'on'); 
demo_hist   = create_fig('fig_A4');  hold(demo_hist, 'on'); 

scale_distr = 0.05; 
x_sclbr = 5.9; y_sclbr = 0.035; 
sclbar_style = {'-k', 'linewidth',0.5};
plot(demo_distr, [1,1]*x_sclbr, [0,scale_distr] + [1,1]*y_sclbr, sclbar_style{:}); 
plot(demo_distr, [0,-0.02] + [1,1]*x_sclbr, [1,1]*y_sclbr, sclbar_style{:}); % to make it pretty 
plot(demo_distr, [0,-0.02] + [1,1]*x_sclbr, [1,1]*(y_sclbr+scale_distr), sclbar_style{:}); % to make it pretty 

raster_ylim = [0, num_unit_trial] + 100*[-1,1]; 
entropy_demo = zeros([1,n_gelec_demo]);
g_elec_ttl = right_lbl_names.g_elec_int_to_int; 
cmap = parula(n_gelec_demo)*0.70;

title_style = {'fontsize', 10, 'fontweight', 'normal'}; 
for i = 1:n_gelec_demo
    loc_i = loc_gelec_demo(i);
    combo_i = {loc_sigma_demo, loc_i};
    cospk_i =  [raw_co_spike{:,combo_i{:}}];
    src_col = cospk_i(1,:);
    int_col = cospk_i(2,:);
    tgt_col = cospk_i(3,:);
   
    [~,order_tgt_src] = sort(tgt_col-src_col);
    if i == 1
        int_uncpld = int_col(1) - src_col(1);
        tgt_uncpld = tgt_col(1) - src_col(1);
    end
    %% Fig 8A1,2
    idx_demo_raster = find(loc_gelec_raster == loc_i); 
    if ~isempty(idx_demo_raster)
        plot(demo_raster(idx_demo_raster), [0,0],raster_ylim, '--', ...
            'linewidth', 2, 'color', cellpop_lbl.Src.color);
        plot(demo_raster(idx_demo_raster), [1,1]*int_uncpld,raster_ylim, ':', ...
            'linewidth', 2, 'color', cellpop_lbl.Int.color);
        plot(demo_raster(idx_demo_raster), [1,1]*tgt_uncpld,raster_ylim, ':', ...
            'linewidth', 2, 'color', cellpop_lbl.Tgt.color);
        for c_i = loc_cell_demo 
            key = cell_keys{c_i};
            dat_col = cospk_i(cellpop_lbl.(key).order,:);
            color_ci = cellpop_lbl.(key).color;
            
            scatter(demo_raster(idx_demo_raster), ...
                dat_col(order_tgt_src) - src_col(order_tgt_src), raster_unit_trial, ...
                marker_size, 'MarkerFaceColor', color_ci, ...
                common_scatter_style{:});            
        end

        xlim(demo_raster(idx_demo_raster), [-1, 25]);
        ylim(demo_raster(idx_demo_raster), raster_ylim);
        set(demo_raster(idx_demo_raster), 'ytick', [1,num_unit_trial], 'box', 'on');             
        ylabel(demo_raster(idx_demo_raster), '# unit \times # trial');
        text(demo_raster(idx_demo_raster), 15, 600,...
            [g_elec_ttl.name ' = ' num2str(g_elec(loc_i)) ' ' g_elec_ttl.unit], ...
            'color','k', 'fontsize', 10, 'fontweight', 'normal');
    end
    
    %% Fig 8A3
    hist_dat = histcounts(tgt_col-src_col, edges);
    hist_dat = hist_dat/num_unit_trial;
    smoothed_hist = conv(hist_dat, smooth_win, 'same');  
    
    fill(demo_distr, centers, smoothed_hist +  0.07*(i-0.5), ...
        cmap(i,:), 'EdgeColor', cmap(i,:), 'LineWidth', 1, 'LineStyle', 'none');
    
    %% Fig 8A4
    ent = -hist_dat.*log2(hist_dat);
    ent(~isfinite(ent)) = 0;
    entropy_demo(i) = sum(ent(:));
    bar(demo_hist, i, entropy_demo(i), 'EdgeColor', 'none', 'FaceColor', cmap(i,:));

    
end

% Additional style for Fig 8A1,2
set(demo_raster(1), 'xtick', '');
xtick_demo = xticks(demo_raster(2)); 
loc_0ms = find(xtick_demo == 0,1); 
xtick_labels = num2cell(xtick_demo); 
xtick_labels{loc_0ms} = ''; 
xticklabels(demo_raster(2),xtick_labels); 

xlabel(demo_raster(2), 'Time (ms)');

% Additional style for Fig 8A3
title(demo_distr, 'Latency (Tgt - Src)', title_style{:}); 
xlabel(demo_distr, 'Time (ms)');
ylabel(demo_distr, 'Density');
set(demo_distr, 'ytick', '', 'xtick',[4.2,6]);
ylim(demo_distr, [0.01, 0.4]);
xlim(demo_distr,[4.2,6]);

% Additional style for Fig 8A4
title(demo_hist, {'Entropy of the Latency (Tgt - Src)'}, title_style{:});
xlabel(demo_hist, '\Sigma G_{elec} (nS)' );
ylabel(demo_hist, 'bits');
set(demo_hist, 'xtick', 1:n_gelec_demo, 'xticklabel',  g_elec(loc_gelec_demo));
ylim(demo_hist, [0, 6]);
%% Fig 8B&C

field2plt = {'mi_src_out', 'mi_over_hout'};  
ann_field2plt = {'B', 'C'}; 

idx_x = 3; 
idx_y = 2; 
x_values = stat_dist{1}.value_var{idx_x}; 
y_values = stat_dist{1}.value_var{idx_y}; 

x_props = right_lbl_names.(stat_dist{1}.var_names{idx_x}); 
y_props = right_lbl_names.(stat_dist{1}.var_names{idx_y}); 

x_range = [min(x_values), max(x_values)];
y_range = [min(y_values), max(y_values)];
dx = abs(x_values(2)-x_values(1)); 
dy = abs(y_values(2)-y_values(1)); 

config_lim = @(v,dv) v + dv * 0.5 * [-1,1];  

loc_pref2plt = strcmp(cell_pref, 'Tgt');

cmap = hot(10000)*0.99; 
colormap(cmap);

for i = 1:length(field2plt)
    ann_i = ann_field2plt{i}; 
    field_i = field2plt{i};
    name_i = right_lbl_names.(field_i).name;
    unit_i = right_lbl_names.(field_i).unit;
    lim_i = right_lbl_names.(field_i).limit;     
    
    create_ann(['ann_' ann_i], ann_i); 
    create_text(['text_' ann_i], name_i, {});
    
    for k = 1:num_cases
        tag_ik = [ann_i, num2str(k)]; 
        ax_plt = create_fig(['fig_' tag_ik]);
        hold(ax_plt, 'on');
        loc2plt = cell(1,3);
        loc2plt{end} = loc_pref2plt;       
        loc2plt{idx_y-1} = ':'; 
        loc2plt{idx_x-1} = ':';
        
        inf_k = stat_dist{k}.info_measures;     
        pix_dat = inf_k(:,:,loc_pref2plt);
        pix_dat = reshape([pix_dat.(field_i)], size(pix_dat));        
        if contains(field_i, 'over')         
            pix_dat = 100*pix_dat;
        end
        
        image(ax_plt, x_values, y_values, pix_dat, 'CDataMapping','scaled');
        caxis(ax_plt, lim_i);         

        xlim(ax_plt, config_lim(x_range, dx)); 
        ylim(ax_plt, config_lim(y_range, dy));

        set(ax_plt, 'fontsize', 10, 'box', 'on', 'linewidth', 1.5); 
        pbaspect(ax_plt, [1,1,1]); 
        
        if k == num_cases            
            clrbar = colorbar(ax_plt); 
            pos_clrbar = layout_map(['colorbar_' ann_i]).normz_pos;
            set(clrbar, 'box', 'off', 'position', pos_clrbar);
            title(clrbar, unit_i, 'fontsize', 10);            
        end
        if k == 1
            ylabel(ax_plt, [y_props.name ' (' y_props.unit ')'], 'fontsize', 12);
        else
            set(ax_plt, 'ytick', '');
        end
        if i == length(field2plt)
            xlabel(ax_plt, [x_props.name ' (' x_props.unit ')'], 'fontsize', 12);
        else
            set(ax_plt, 'xtick', '');
        end
    end
end
