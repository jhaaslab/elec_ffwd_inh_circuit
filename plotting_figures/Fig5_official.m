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
    
    data_file_pref = dir(fullfile(data_folder, filename_prefix));
    data_file_pref = {data_file_pref.name};
    data2save = struct();
    for i = 1:length(data_file_pref)
        loaded_data = load(fullfile(data_folder,data_file_pref{i}));
        fN = fieldnames(loaded_data.data2save);
        for j = 1:length(fN)
            data2save.(fN{j}) = loaded_data.data2save.(fN{j});
        end
    end
    
    vars_ = data2save.vars;
    var_names = fieldnames(vars_);
    n_vars = length(var_names);
    order_var = cellfun(@(x) double(vars_.(x).order), var_names,'UniformOutput',true);
    unit_var = cellfun(@(x) vars_.(x).unit, var_names,'UniformOutput',false);
    value_var = cellfun(@(x) double(vars_.(x).range), var_names,'UniformOutput',false);
    len_var = cellfun(@(x) length(x), value_var, 'UniformOutput', true);
    
    cell_names = cellstr(data2save.cell_names)';
    cell_pref = {'Src', 'Int', 'Tgt'};
    num_pref = length(cell_pref);
    idx_pref  = struct();
    for i = 1:num_pref
        idx_pref.(cell_pref{i}) = find(contains(cell_names, cell_pref{i}));
    end
    
    raw_data = cell([len_var',num_pref]);
    t_inp = cell(len_var');
    for i = 1:data2save.num_dat
        dat_ = data2save.(sprintf('dat_%06d',i));
        var_ = double(dat_.var);
        loc_ = cell([1, n_vars]);
        for k = 1:n_vars
            loc_k = find(abs(value_var{k} - var_(order_var(k)))<eps);
            if isempty(loc_k)
                error('Could not find ''index'' of %s at location %d', ...
                    sprintf('dat_%06d',i), k);
            end
            loc_{k} = loc_k;
        end
        index = double(dat_.dat.idx_v)+1;
        spkts = dat_.dat.t_v;
        t_inp{loc_{:}} = dat_.ti_inp;
        for j = 1:num_pref
            index2find = idx_pref.(cell_pref{j});
            spkt_j = [];
            for idx = index2find
                spkt_j = [spkt_j, spkts(index == idx)];  %#ok<AGROW>
            end
            raw_data{loc_{:},j} = spkt_j;
        end
    end
    
    %% Normalize data
    sizeRD = size(raw_data);
    spkt_dist = cell(sizeRD(2:end));
    tinp_dist = cell(sizeRD(2:end-1));
    stat_dim = num2cell(sizeRD(2:end));
    clear stat_dist;
    stat_dist(stat_dim{:}) = struct('mean', [], 'std', [], ...
        'ptsh', [], 'centers', [], ...
        'peak', [], 'auc', [], ...
        'lat', [], 'lat5per', [], 'lat10per', []);
    
    
    grid_arrays = arrayfun(@(x) {1:x}, sizeRD(2:end-1));
    grid_cells = cell(1,length(grid_arrays));
    [grid_cells{:}] = ndgrid(grid_arrays{:});
    tot_combos = numel(grid_cells{1}(:));
    grid_combos = zeros(tot_combos, length(grid_arrays));
    for i = 1:length(grid_arrays)
        grid_combos(:,i) = grid_cells{i}(:);
    end
    Src_loc = find(strcmp(cell_pref, 'Src'));
    f_binw = @(x) (x/10);
    t_init_std = value_var{2};
    edge_lim = [0, 200];
    smw = 20;
    n_trials = sizeRD(1);
    
    for i = 1:tot_combos
        combo_i = num2cell(grid_combos(i,:));
        
        binw = f_binw(t_init_std(combo_i{1}));
        n_bins = ceil(abs(diff(edge_lim))/binw);
        Src_dat = [raw_data{:,combo_i{:},Src_loc}];
        [ptsh_src, centers_src] = return_histogram(Src_dat, n_trials, n_bins, edge_lim, smw);
        mean_src = mean(Src_dat);
        std_src = std(Src_dat);
        peak_ptsh_src = max(ptsh_src);
        auc_ptsh_src = trapz(centers_src, ptsh_src);
        
        for j = 1:sizeRD(end)
            dat_combo_ij = [raw_data{:,combo_i{:},j}];
            spkt_dist{combo_i{:},j} = [raw_data{:,combo_i{:},j}];
            
            stat_dist(combo_i{:},j).mean = mean(dat_combo_ij)-mean_src;
            stat_dist(combo_i{:},j).std = std(dat_combo_ij)/std_src;
            
            [ptsh_comboij, centers_comboij] = return_histogram(dat_combo_ij, n_trials, n_bins, edge_lim, smw);
            normz_ptsh = ptsh_comboij/peak_ptsh_src;
            normz_centers = centers_comboij-mean_src;
            peak_normz_ptsh = max(normz_ptsh);
            t_peak = normz_centers(find(normz_ptsh == peak_normz_ptsh,1));
            t_5percent_srcpeak = max(normz_centers(normz_ptsh < 0.05 & normz_centers < t_peak));
            t_10percent_srcpeak = max(normz_centers(normz_ptsh < 0.10 & normz_centers < t_peak));
            
            stat_dist(combo_i{:},j).ptsh = normz_ptsh;
            stat_dist(combo_i{:},j).centers = normz_centers;
            stat_dist(combo_i{:},j).peak = peak_normz_ptsh;
            stat_dist(combo_i{:},j).auc = trapz(centers_comboij, ptsh_comboij)/auc_ptsh_src;
            
            %% all latencies are normalized to src_stddev
            stat_dist(combo_i{:},j).lat = normz_centers(find(normz_ptsh>0,1))/std_src;
            stat_dist(combo_i{:},j).lat5per = t_5percent_srcpeak/std_src;
            stat_dist(combo_i{:},j).lat10per = t_10percent_srcpeak/std_src; %% we used this 
        end
        tinp_dist{combo_i{:}} = [t_inp{:,combo_i{:}}];
    end
end

%% Load layout
file_name  = 'Fig5_layout_v2.svg';
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

%% Value vectors
t_init_std = value_var{2};
g_elec_vec = value_var{3};

f_vec2show = @(x) [x(1), x(ceil(length(x)/2)), x(end)];
name_x_ = '\sigma_{inp}'; 
x_axis_actual = struct(...
    'name', [name_x_ '(ms)'], ...
    'vec', t_init_std, ...
    'show', f_vec2show(t_init_std) );
y_axis_actual = struct(...
    'name', '\SigmaG_{elec} (nS)', ...
    'vec', g_elec_vec, ...
    'show', f_vec2show(g_elec_vec) );

%% Fig 5A
create_ann('ann_A', 'A');

%% Text for Fig 5B
add_text_style = {'fontweight', 'normal', 'fontsize', 16}; 
create_text('text_B@1', 'Int', add_text_style); 
create_text('text_B@2', 'Tgt', add_text_style); 

%% Fig 5B
f_limx = @(x) [-0.5,1]*4*x + [0,5];
f_limx_inset = @(x) [-1,-1]*3*x + [+4,5+x/1.12];
f_limy = @(y) [-0.1, 2.5];
create_ann('ann_B', 'B');


std_tinit2plt = [1, 3, 5, 10];
cells2plt = {'Int', 'Tgt'};

nline_plt = length(g_elec_vec);
cmap = parula(nline_plt);
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
                'color', cmap(k,:), 'linewidth', 1);
            src_ptsh = [src_ptsh; stat_dist(loc_i,k,loc_src).ptsh]; %#ok<AGROW>
        end
        if j == 1
            plot(ax_ij, stat_dist(loc_i,1,loc_src).centers,...
                mean(src_ptsh,1), '--k', 'linewidth', 1);
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
            set(ax_ij, 'ytick', ''); 
        end 
        set(ax_ij, 'TickLength', [0.025,0.025]); 
        
        
        
        %% inset to zoom in Tgt latency
        if strcmp(cells2plt{j}, 'Tgt')
            ax_inset = create_fig(['inset_B' num2str(i) '_' num2str(j)]) ;
            hold(ax_inset, 'on');
            for k = 1:nline_plt
                stat_dist_ikj = stat_dist(loc_i,k,loc_j);
                plot(ax_inset, stat_dist_ikj.centers, stat_dist_ikj.ptsh, ...
                    'color', cmap(k,:), 'linewidth', 0.5);
            end
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
