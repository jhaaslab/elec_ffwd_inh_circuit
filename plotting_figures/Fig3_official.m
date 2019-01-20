clc; clear; close all;
addpath(genpath('../figures'), genpath('../functions'));
%% Load data
data_folder = '../data/net2_sim/';
data_file_pref = dir(fullfile(data_folder,'Net2_GJvsdtInp_0*.mat'));
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

Src1_spkt = [];
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
    
    if isempty(Src1_spkt)
        [~,loc_spk] = findpeaks(dat_.dat.Src_1, t, 'MinPeakProminence', 50);
        Src1_spkt = loc_spk;
    end
    
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
file_name  = 'Fig3_layout_official.svg';
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
%% Fig 3A
create_ann('ann_A', 'A');

%% Fig 3B
create_text('text_B', 'PSP of Tgt', {'FontSize', 15}); 
tscale_start = 70.5;
pspscale_start = -2;
time_scale = 1;
psp_scale = 1;
 
int_line_clr = 0.4*[1,1,1]; 

dt2find = [0.6,1.2,2.4,4.8]; 
x_loc2plt = arrayfun(@(rr) find(abs(x_vec-rr)<eps*1e3),dt2find,'UniformOutput',true); 
y_loc2plt = 1:15:length(y_vec);
cmap = parula(length(y_loc2plt)) * gelec_cmap_factor;
nsplt = ceil(sqrt(length(x_loc2plt)));

for idx_x = 1:length(x_loc2plt)
    tag_idx = num2str(idx_x);
    if idx_x == 1 
        create_ann(['ann_B' tag_idx], 'B');
    end

    ax_spkt = create_fig(['fig_B' tag_idx '_2']);
    hold(ax_spkt, 'on');
        
    ax_tgt = create_fig(['fig_B' tag_idx '_1']);
    hold(ax_tgt, 'on');
    plot(ax_tgt, [0 100], [0, 0], ':', 'color', 0.6*[1,1,1], 'linewidth', 0.75);
    set(ax_tgt, 'visible', 'off') 
    
    ax_int_line =  create_fig(['fig_B' tag_idx '_3']);
    hold(ax_int_line, 'on'); 
    
    i_x = x_loc2plt(idx_x);
    
    plot(ax_int_line, y_vec, [Int1_spkt{i_x,:}] - Src1_spkt, '-', 'linewidth', 1, 'Color', int_line_clr);     
    plot(ax_int_line, y_vec, [Int2_spkt{i_x,:}] - Src1_spkt, ':', 'linewidth', 1.5,  'Color', int_line_clr);      
    x_lim = [floor(y_vec(1)), ceil(y_vec(end))]; 
    y_lim = [0, 8];  
    xlim(ax_int_line, x_lim); 
    ylim(ax_int_line, y_lim);    
    set(ax_int_line, 'xtick', x_lim, 'ytick', y_lim , ...
        'fontsize', 8, 'color', 'none', ...
        'xcolor', int_line_clr, 'ycolor', int_line_clr); 
    if (idx_x == 1) 
        xlabel(ax_int_line, 'G_{elec} (nS)'); 
        ylabel(ax_int_line, 'Latency (ms)');
    end
        
    for idx_y = 1:length(y_loc2plt)
        i_y = y_loc2plt(idx_y);
        tgt_xy = reshape(Tgt_PSP(:,i_x,i_y), [lenT, 1]);
        int1_xy = Int1_spkt{i_x,i_y};
        int2_xy = Int2_spkt{i_x,i_y};
        clr_xy = cmap(idx_y,:);
        
        plot(ax_tgt, t, tgt_xy, '-', 'Color', clr_xy, 'LineWidth', 1);
        
        ycoord_spkt = -0.75 + (1.35+0.75)*idx_y/length(y_loc2plt);
        spkt_marker_prop = {'MarkerSize', 5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', clr_xy}; 
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
    
    xlim(ax_tgt, [70, 80]);
    ylim(ax_tgt, [-3, 5]);
    plot(ax_tgt, tscale_start+[0,time_scale], pspscale_start*[1,1], ...
        '-k', 'LineWidth', 0.75);
    plot(ax_tgt, tscale_start*[1,1], pspscale_start+[0,psp_scale], ...
        '-k', 'LineWidth', 0.75);
    title(ax_tgt, ['\Deltat_{inp} = ' num2str(x_vec(i_x)) ' ms'],...
        'fontweight', 'normal', 'fontsize', 10, 'visible', 'on', 'HorizontalAlignment', 'right'); 
    set(ax_tgt, 'xcolor', 'none', 'ycolor', 'none');
    
    xlim(ax_spkt, [70 80]);
    ylim(ax_spkt, [-1.5, 1.45]);
    set(ax_spkt, 'visible', 'off');
end

colormap(ax_tgt, cmap(1:end,:));
clrbar = colorbar(ax_tgt);
colorbar_show = y_vec(y_loc2plt);
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', linspace(0,1,length(colorbar_show)), 'ticklabels', colorbar_show, ...
    'Position', pos_clrbar, 'fontsize', 9);
xlabel(clrbar, 'G_{elec} (nS)', 'fontsize', 10);

%% Fig 3 CDE
splt_tag = {'C', 'D', 'E'};
field2plt = {'peak', 'iw', 'auc'};

axislabel_style = {'fontsize',8}; 
[xvm, yvm] = meshgrid(x_vec, y_vec); 
for i = 1:length(splt_tag)
    tag_i = splt_tag{i};
    field_i = psp_anly.(field2plt{i});
    data2plt_i = field_i.dat; 
    title_i = field_i.name; 
    unit_i = field_i.unit; 
    create_ann(['ann_' tag_i '1'], [tag_i '_1']);
    create_ann(['ann_' tag_i '2'], [tag_i '_2']);
    
    %% CDE_2: pixels 
    ax_pixels = create_fig(['fig_' tag_i '2']);
    hold(ax_pixels, 'on');
    image(ax_pixels,x_vec, y_vec, data2plt_i', 'CDataMapping','scaled');
    set(ax_pixels, 'ydir', 'normal')
    xlabel(ax_pixels, '\Deltat_{inp} (ms)', axislabel_style{:}); 
    ylabel(ax_pixels, 'G_{elec} (nS)', axislabel_style{:}); 
    xlim(ax_pixels, [min(x_vec),max(x_vec)]);
    ylim(ax_pixels, [min(y_vec),max(y_vec)]);
    set(ax_pixels, 'box', 'on', 'linewidth', 1.5); 
    colormap(ax_pixels, 'hot'); 
    cmap = colormap(ax_pixels); 
    cmap(end,:) = ''; 
    colormap(ax_pixels, cmap);
    clrbar = colorbar(ax_pixels); 
    caxis(ax_pixels, [min(data2plt_i(:)),max(data2plt_i(:))]);
    title(clrbar, unit_i); 
    set(clrbar, 'Box', 'off');

    %% CDE_1: lines 
    create_text(['text_' tag_i '1'], title_i, {}); 
    ax_lines = create_fig(['fig_' tag_i '1']);
    hold(ax_lines, 'on');
    cmap = parula(length(y_vec)) * gelec_cmap_factor;
    for j = 1:length(y_vec)
        plot(ax_lines, x_vec, data2plt_i(:,j), '-', 'LineWidth', 0.5, 'Color', cmap(j,:));
    end
    xlabel(ax_lines, '\Deltat_{inp} (ms)', axislabel_style{:}); 
    ylabel(ax_lines, unit_i, axislabel_style{:}); 
    xlim(ax_lines, [min(x_vec),max(x_vec)]);
    colormap(ax_lines, cmap); 
    clrbar = colorbar(ax_lines); 
    set(clrbar, 'Box', 'off')
    title(clrbar, 'G_{elec}'); 
    caxis(ax_lines, [min(y_vec), max(y_vec)]); 
end