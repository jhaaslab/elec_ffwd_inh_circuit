clc; clear; close all;
addpath(genpath('../figures'), genpath('../functions'));
%% Load data
data_folder = '../data/net2_sim/';
data_file_pref = dir(fullfile(data_folder,'Net2_GJvsdtInp_Normz*.mat'));
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

load('../data/net1_sim/ctrl_anly.mat');
vars_ = data2save.vars;
var_names = fieldnames(vars_);
order_var = cellfun(@(x) double(vars_.(x).order), var_names,'UniformOutput',true);
unit_var = cellfun(@(x) vars_.(x).unit, var_names,'UniformOutput',false);
value_var = cellfun(@(x) vars_.(x).range, var_names,'UniformOutput',false);
x_name = var_names{order_var(2)};
x_vec = value_var{order_var(2)};
x_unit = unit_var{order_var(2)};
y_name = var_names{order_var(1)};
y_vec = value_var{order_var(1)};
y_unit = unit_var{order_var(1)};

t = data2save.dat_000001.dat.t;
lenT = length(t);

Int1_spkt = cell(length(x_vec),length(y_vec));
Int2_spkt = cell(length(x_vec),length(y_vec));
Tgt_PSP = zeros(lenT,length(x_vec),length(y_vec));
iw_psp = zeros(length(x_vec),length(y_vec));
peak_psp = zeros(length(x_vec),length(y_vec));
auc_psp = zeros(length(x_vec),length(y_vec));

for i = 1:data2save.num_dat
    dat_ = data2save.(sprintf('dat_%06d',i));
    var_ = dat_.var;
    x_loc = find(abs(x_vec-var_(order_var(1)))<eps);
    y_loc = find(abs(y_vec-var_(order_var(2)))<eps);
    if isempty(x_loc) && isempty(y_loc)
        error('Could not find ''index'' of %s', sprintf('dat_%06d',i));
    end

    tgt_vm = dat_.dat.Tgt;
    Tgt_PSP(:,x_loc,y_loc) = tgt_vm - tgt_vm(1);

    [~,loc_spk] = findpeaks(dat_.dat.Int_1, t, 'MinPeakProminence', 50);
    Int1_spkt{x_loc,y_loc} = loc_spk;

    [~,loc_spk] = findpeaks(dat_.dat.Int_2, t, 'MinPeakProminence', 50);
    Int2_spkt{x_loc,y_loc} = loc_spk;

    iw_psp(x_loc,y_loc) = dat_.iw.dt_epsp;
    peak_psp(x_loc,y_loc) = dat_.iw.peak_epsp;
    auc_psp(x_loc,y_loc) = dat_.iw.auc_epsp;
end
psp_anly = struct(  'iw', struct('dat', iw_psp, 'name', 'Integration window', 'unit', 'ms'), ...
                    'peak', struct('dat', peak_psp, 'name', 'Peak of PSP', 'unit', 'mV'), ...
                    'auc', struct('dat', auc_psp, 'name', 'AUC of PSP', 'unit', 'ms\cdotmV'));
clear iw_psp peak_psp auc_psp;

%% Load layout
file_name  = 'Fig3S_layout.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor);

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

%% Fig 3A
create_ann('ann_A', 'A');

%% Fig 3B
create_text('text_B', 'PSP of Tgt', {}); 
tscale_start = 70.5;
pspscale_start = -2;
time_scale = 1;
psp_scale = 1;
ax_spkt_bckgrnd = 0.9*[1,1,1];
dt2find = [1.5, 3, 4.5, 8, 12, 18]; 
x_loc2plt = arrayfun(@(rr) find(abs(x_vec-rr)<eps*1e3),dt2find,'UniformOutput',true); 
y_loc2plt = 1:15:length(y_vec);
clr_fctr = 1;
cmap = parula(clr_fctr*length(y_loc2plt));
nsplt = ceil(sqrt(length(x_loc2plt)));

ctrl_t = ctrl_anly.t;
for idx_x = 1:length(x_loc2plt)
    tag_idx = num2str(idx_x);
    if idx_x == 1 
        create_ann(['ann_B' tag_idx], 'B');
    end
    ax_tgt = create_fig(['fig_B' tag_idx '_1']);
    hold(ax_tgt, 'on');
    plot(ax_tgt, [0 100], [0, 0], ':', 'color', 0.6*[1,1,1], 'linewidth', 0.75);
    ax_spkt = create_fig(['fig_B' tag_idx '_2']);
    hold(ax_spkt, 'on');
    
    i_x = x_loc2plt(idx_x);
    for idx_y = 1:length(y_loc2plt)
        i_y = y_loc2plt(idx_y);
        tgt_xy = reshape(Tgt_PSP(:,i_x,i_y), [lenT, 1]);
        int1_xy = Int1_spkt{i_x,i_y};
        int2_xy = Int2_spkt{i_x,i_y};
        clr_xy = cmap(clr_fctr*idx_y,:);
        
        plot(ax_tgt, t, tgt_xy, '-', 'Color', clr_xy, 'LineWidth', 1);
        
        ycoord_spkt = -0.75 + (1.35+0.75)*idx_y/length(y_loc2plt);
        spkt_marker_prop = {'MarkerSize', 4, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', clr_xy}; 
        for i1 = 1:length(int1_xy)
            y_i1 = ycoord_spkt; 
            plot(ax_spkt, int1_xy(i1), y_i1, ...
                'Marker', 'o', spkt_marker_prop{:});
        end
        for i2 = 1:length(int2_xy)
            y_i2 = ycoord_spkt; 
            plot(ax_spkt, int2_xy(i2), y_i2, ...
                'Marker', '^', spkt_marker_prop{:});
        end
    end
    
    loc_ctrl = abs(ctrl_anly.dt_inp - x_vec(i_x))<eps;
    ctrl_tgt = ctrl_anly.Tgt(:,loc_ctrl);
    ctrl_tgt = ctrl_tgt - ctrl_tgt(1);
    plot(ax_tgt, ctrl_t, ctrl_tgt, '-.k');
    [~,int_ctrl] = findpeaks(ctrl_anly.Int(:,loc_ctrl), ctrl_t, 'MinPeakProminence', 50);
    for i_c = 1:length(int_ctrl)
        plot(ax_spkt, int_ctrl(i_c), -1, 'ok', 'MarkerSize',3.5)
    end
    
    xlim(ax_tgt, [70 100]);
    ylim(ax_tgt, [-5, 4]);
    plot(ax_tgt, tscale_start+[0,time_scale], pspscale_start*[1,1], ...
        '-k', 'LineWidth', 0.75);
    plot(ax_tgt, tscale_start*[1,1], pspscale_start+[0,psp_scale], ...
        '-k', 'LineWidth', 0.75);
    title(ax_tgt, ['\Deltat_{inp} = ' num2str(x_vec(i_x)) 'ms'], 'fontweight', 'normal')
    set(ax_tgt, 'xcolor', 'none', 'ycolor', 'none');
    
    xlim(ax_spkt, [70 100]);
    ylim(ax_spkt, [-1.5, 1.45]);
    set(ax_spkt, 'visible', 'off');
end

colormap(ax_tgt, cmap(clr_fctr:clr_fctr:end,:));
clrbar = colorbar(ax_tgt);
colorbar_show = y_vec(y_loc2plt(1:end));
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', linspace(0,1,length(colorbar_show)), 'ticklabels', colorbar_show, ...
    'Position', pos_clrbar);
title(clrbar, 'G_{elec} (nS)');

%% Fig 3 CDE
splt_tag = {'C', 'D', 'E'};
field2plt = {'peak', 'iw', 'auc'};
corr_ctrl_attr = {'peak_epsp', 'dt_epsp', 'auc_epsp'};

axislabel_style = {'fontsize',7}; 
[xvm, yvm] = meshgrid(x_vec, y_vec);
ctrl_dtinp = ctrl_anly.dt_inp; 

for i = 1:length(splt_tag)
    tag_i = splt_tag{i};
    field_i = psp_anly.(field2plt{i});
    data2plt_i = field_i.dat; 

    title_i = field_i.name; 
    unit_i = field_i.unit; 
    ctrl_dat_i = ctrl_anly.(corr_ctrl_attr{i}); 
    create_ann(['ann_' tag_i '1'], [tag_i '_1']);
    create_ann(['ann_' tag_i '2'], [tag_i '_2']);
    
    %% CDE_1 
    ax_pixels = create_fig(['fig_' tag_i '2']); 
    hold(ax_pixels, 'on');
    image(ax_pixels,x_vec, y_vec, data2plt_i', 'CDataMapping','scaled');
    set(ax_pixels, 'ydir', 'normal')
    xlabel(ax_pixels, '\Deltat_{inp} (ms)', axislabel_style{:}); 
    ylabel(ax_pixels, 'G_{elec} (nS)', axislabel_style{:}); 
    xlim(ax_pixels, [min(x_vec),max(x_vec)]);
    ylim(ax_pixels, [min(y_vec),max(y_vec)]);
    
    colormap(ax_pixels, 'hot'); 
    cmap = colormap(ax_pixels); 
    cmap(end,:) = ''; 
    colormap(ax_pixels, cmap);
    clrbar = colorbar(ax_pixels); 
    caxis(ax_pixels, [min(data2plt_i(:)),max(data2plt_i(:))]);
    title(clrbar, unit_i); 
    set(clrbar, 'Box', 'off', 'LineWidth', 0.001)

    %% CDE_2
    create_text(['text_' tag_i '1'], title_i, {});
    ax_lines = create_fig(['fig_' tag_i '1']); 
    hold(ax_lines, 'on');
    y_vec2plt = y_vec(1:5:end); 
    data2plt_i = data2plt_i(:,1:5:end); 
    cmap = parula(length(y_vec2plt));
    for j = 1:length(y_vec2plt)
        plot(ax_lines, x_vec, data2plt_i(:,j), '-', 'LineWidth', 0.5, 'Color', cmap(j,:));
    end
    plt_ctrl = plot(ctrl_dtinp, ctrl_dat_i, '-.k', 'LineWidth', 1.5); 
    xlabel(ax_lines, '\Deltat_{inp} (ms)', axislabel_style{:}); 
    ylabel(ax_lines, unit_i, axislabel_style{:}); 
    xlim(ax_lines, [min(x_vec),max(x_vec)]);
    ylim(ax_lines, [min(data2plt_i(:)) - 0.25, max(data2plt_i(:)) + 0.25]); 
    colormap(ax_lines, cmap); 
    clrbar = colorbar(ax_lines); 
    set(clrbar, 'Box', 'off', 'LineWidth', 0.001)
    title(clrbar, 'G_{elec}'); 
    caxis(ax_lines, [min(y_vec), max(y_vec)]); 
    
end