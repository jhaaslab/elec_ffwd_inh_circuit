clc; clear; close all;
addpath(genpath('../figures'), genpath('../functions'));

%% Load data
data_folder = '../data/net3_sim/';
data_file_pref = dir(fullfile(data_folder,'Net3_GJvsGABAIntvsdtInp_*.mat'));
data_file_pref = {data_file_pref.name};
data2save = struct();
for i = 1:length(data_file_pref)
    loaded_data = load(fullfile(data_folder,data_file_pref{i}));
    fN = fieldnames(loaded_data.data2save);
    for j = 1:length(fN)
        data2save.(fN{j}) = loaded_data.data2save.(fN{j});
    end
end

for i = 1:data2save.num_dat
    nme_ = sprintf('dat_%06d',i);
    dat_ = data2save.(nme_);
    data2save.(nme_).iw = analyze_TgtPSP(dat_.dat.t, dat_.dat.Tgt);
end

vars_ = data2save.vars;
var_names = fieldnames(vars_);
n_vars = length(var_names);
order_var = cellfun(@(x) double(vars_.(x).order), var_names,'UniformOutput',true);
unit_var = cellfun(@(x) vars_.(x).unit, var_names,'UniformOutput',false);
value_var = cellfun(@(x) double(vars_.(x).range), var_names,'UniformOutput',false);
len_var = cellfun(@(x) length(x), value_var, 'UniformOutput', true);

t = data2save.dat_000001.dat.t;
lenT = length(t);

Int1_spkt = cell(len_var');
Int2_spkt = cell(len_var');
Tgt_PSP = zeros([lenT,len_var']);
iw_psp = zeros(len_var');
peak_psp = zeros(len_var');
auc_psp = zeros(len_var');


for i = 1:data2save.num_dat
    dat_ = data2save.(sprintf('dat_%06d',i));
    var_ = dat_.var;
    loc_ = cell([1, n_vars]);
    for k = 1:n_vars
        loc_k = find(abs(value_var{k} - var_(order_var(k)))<eps);
        if isempty(loc_k)
            error('Could not find ''index'' of %s at %d location', ...
                sprintf('dat_%06d',i), k);
        end
        loc_{k} = loc_k;
    end

    tgt_vm = dat_.dat.Tgt;
    Tgt_PSP(:,loc_{:}) = tgt_vm - tgt_vm(1);

    [~,loc_spk] = findpeaks(dat_.dat.Int_1, t, 'MinPeakProminence', 50);
    Int1_spkt{loc_{:}} = loc_spk;

    [~,loc_spk] = findpeaks(dat_.dat.Int_2, t, 'MinPeakProminence', 50);
    Int2_spkt{loc_{:}} = loc_spk;

    iw_psp(loc_{:}) = dat_.iw.dt_epsp;
    peak_psp(loc_{:}) = dat_.iw.peak_epsp;
    auc_psp(loc_{:}) = dat_.iw.auc_epsp;
end
psp_anly = struct(  'iw', struct('dat', iw_psp, 'name', 'Integration window', 'unit', 'ms'), ...
                    'peak', struct('dat', peak_psp, 'name', 'Peak of PSP', 'unit', 'mV'), ...
                    'auc', struct('dat', auc_psp, 'name', 'AUC of PSP', 'unit', 'ms\cdotmV'));
clear iw_psp peak_psp auc_psp;

%% Load layout
file_name  = 'Fig4_layout_official.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor, ...
    'PaperUnits', 'centimeters','PaperPosition', [0, 0, width, height]*conv_factor,...
    'PaperSize', [width, height]*conv_factor);
%% Annotation and axes styles
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 15, 'FontWeight', 'bold'};
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

gelec_cmap_factor = 0.85; 
%% Fig 4A
create_ann('ann_A', 'A');

%% Value vectors
dt_inp = value_var{1};
g_gabaint = value_var{2};
g_elecint = value_var{3};

%% Fig 4B
create_text('text_B', 'PSP of Tgt', {}); 
tscale_start = 40.5;
pspscale_start = -2;
time_scale = 1;
psp_scale = 1;

g_gabaint_2look4v = [1, 7];
dt_inp_2look4v = [2, 4.4];
tag_order = [1, 3; 2, 4];
g_elec2show = 1:2:length(g_elecint);
cmap = parula(length(g_elec2show)) * gelec_cmap_factor;

for m = 1:length(dt_inp_2look4v)
    loc_ = cell(1,length(len_var));
    dt_inp_2look4 = dt_inp_2look4v(m);
    loc_{1} = find(abs(dt_inp - dt_inp_2look4) < eps);
    for i = 1:length(g_gabaint_2look4v)
        g_gabaint_2look4 = g_gabaint_2look4v(i);
        
        loc_{2} = find(abs(g_gabaint - g_gabaint_2look4) < eps);
        tag_i = num2str(tag_order(m,i));
        if tag_order(m,i) == 1
            create_ann('ann_B1', 'B');
        end
        ax_spkt = create_fig(['fig_B' tag_i '_2']);
        hold(ax_spkt, 'on');
        
        ax_tgt = create_fig(['fig_B' tag_i '_1']);
        set(ax_tgt, 'visible', 'off'); 
        hold(ax_tgt, 'on');
        plot(ax_tgt, [0 100], [0, 0], ':', 'color', 0.6*[1,1,1], 'linewidth', 0.75);
        for j = 1:length(g_elec2show)
            loc_{3} = g_elec2show(j);
            tgt_ij = reshape(Tgt_PSP(:,loc_{:}), [lenT, 1]);
            int1_ij = Int1_spkt{loc_{:}};
            int2_ij = Int2_spkt{loc_{:}};
            clr_ij = cmap(j,:);
            
            plot(ax_tgt, t, tgt_ij, '-', 'Color', clr_ij, 'LineWidth', 1);
            
            ycoord_spkt = j/length(g_elec2show)-0.1;
            spkt_marker_prop = {'MarkerSize',5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', clr_ij}; 
            for i1 = 1:length(int1_ij)
                y_i1 = ycoord_spkt; 
                plot(ax_spkt, int1_ij(i1), y_i1, ...
                    'Marker', 'o', spkt_marker_prop{:});
            end
            for i2 = 1:length(int2_ij)
                y_i2 = ycoord_spkt; 
                plot(ax_spkt, int2_ij(i2), y_i2, ...
                    'Marker', '^', spkt_marker_prop{:});
            end
        end
        xlim(ax_tgt, [40, 50]);
        ylim(ax_tgt, [-3, 3.5]);
        set(ax_tgt, 'visible', 'off');
        xlim(ax_spkt, [40, 50]);
        ylim(ax_spkt, [-0.2, 1.2]);
        set(ax_spkt, 'visible', 'off');
        
        plot(ax_tgt, tscale_start+[0,time_scale], pspscale_start*[1,1], ...
            '-k', 'LineWidth', 0.75);
        plot(ax_tgt, tscale_start*[1,1], pspscale_start+[0,psp_scale], ...
            '-k', 'LineWidth', 0.75);
        
        title(ax_tgt, ['\Deltat_{inp} = ' num2str(dt_inp_2look4) ' ms'],...
        'fontweight', 'normal', 'fontsize', 10, 'visible', 'on', 'HorizontalAlignment', 'right');
    end
end
colormap(ax_tgt, cmap);
clrbar = colorbar(ax_tgt);
colorbar_show = g_elecint(g_elec2show);
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', linspace(0,1,length(colorbar_show)), 'ticklabels', colorbar_show, ...
    'Position', pos_clrbar);
title(clrbar, 'G_{elec} (nS)');

%% Fig 4CD
splt_tag = {'C', 'D'};
field2plt = {'iw', 'auc'};

g_gabaint_2look4v = [1,3,5,7];

cmap = parula(length(g_elecint)) * gelec_cmap_factor;
for m = 1:length(splt_tag)
    field_m = psp_anly.(field2plt{m});
    data2plt_m = field_m.dat;
    title_m = field_m.name;
    unit_m = field_m.unit;
    lim_y = [min(data2plt_m(:))-1, max(data2plt_m(:))+1];
    create_text(['text_' splt_tag{m}], title_m, {});
    for i = 1:length(g_gabaint_2look4v)
        tag_i = [splt_tag{m} num2str(i)];
        
        g_gabaint_2look4 = g_gabaint_2look4v(i);
        loc_ = cell(1,length(len_var));
        loc_{1} = ':';
        loc_{2} = find(abs(g_gabaint - g_gabaint_2look4) < eps);
        if i == 1
            create_ann(['ann_' splt_tag{m}], splt_tag{m});
        end
        ax_lines = create_fig(['fig_' tag_i]);
        hold(ax_lines, 'on');
        
        for j = 1:length(g_elecint)
            loc_{3} = j;
            data2plt_ij = data2plt_m(loc_{:});
            clr_ij = cmap(j,:);
            
            plot(ax_lines, dt_inp, data2plt_ij, '-', 'Color', clr_ij, 'LineWidth', 1.5);
        end
        xlim(ax_lines, [0, max(dt_inp)]);
        ylim(ax_lines, lim_y);
        xlabel(ax_lines, '\Deltat_{inp} (ms)');
        if i == 1
            ylabel(ax_lines, unit_m);
        end
    end
end
colormap(ax_lines, cmap);
clrbar = colorbar(ax_lines);
colorbar_show = g_elecint(1:2:end);
pos_clrbar = layout_map('colorbar_CD').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', linspace(0,1,length(colorbar_show)), 'ticklabels', colorbar_show, ...
    'Position', pos_clrbar);
title(clrbar, 'G_{elec} (nS)', 'fontsize', 10);
