clc; clear; close all;
addpath(genpath('../figures'), genpath('../functions')); 
%% Load data
data_folder = '../data/net2_sim/';
data_file_pref = dir(fullfile(data_folder,'Net2_GABAvsdtInp_0*.mat'));
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
Src2_spkt = cell(length(x_vec),length(y_vec));
Int1_spkt = [];
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
        
        [~,loc_spk] = findpeaks(dat_.dat.Int_1, t, 'MinPeakProminence', 50);
        Int1_spkt = loc_spk;
    end
    
    [~,loc_spk] = findpeaks(dat_.dat.Src_2, t, 'MinPeakProminence', 50); 
    Src2_spkt{x_loc,y_loc} = loc_spk;
    
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
file_name  = 'Fig2_layout_official.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor, ...
    'PaperUnits', 'centimeters','PaperPosition', [0, 0, width, height]*conv_factor, 'PaperSize', [width, height]*conv_factor);

%% Annotation and axes styless
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 15, 'FontWeight', 'bold'}; 

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:}); 
create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);
title_size = 12; 

%% Fig 2A
create_ann('ann_A', 'A');

%% Fig 2B
create_ann('ann_B', 'B');

time_scale = 1; %ms
psp_scale = 1; %mV 

gaba2plt = 8; 
x_loc2plt = 1:5:length(x_vec); 
y_loc2plt = find(abs(y_vec - gaba2plt) < eps); 

ax_spkt = create_fig('fig_B_2');
hold(ax_spkt, 'on'); 
plot(ax_spkt, Src1_spkt*[1,1], [-2.5,1.2], '-k', 'LineWidth', 1.5); 
text(ax_spkt, Src1_spkt-0.35, -0.45, 'Src_1', 'Fontsize', 10) ;
plot(ax_spkt, Int1_spkt*[1,1], [-2.5,1.2], '-', 'color', 0.6*[1,1,1], 'LineWidth', 1.5); 
text(ax_spkt, Int1_spkt-0.35, -0.45, 'Int_1', 'Fontsize', 10) ;
          
ax_tgt = create_fig('fig_B_1'); 
hold(ax_tgt, 'on'); 
plot(ax_tgt, [0, 200], [0, 0], ':', 'color', 0.6*[1,1,1], 'linewidth',0.75);
gaba_txt = ['G_{GABA\rightarrowTgt} = '  num2str(gaba2plt) 'nS']; 
% text(ax_tgt, 69.5, -1.5, gaba_txt, 'fontsize', 10);

num_lines = length(x_loc2plt); 
cmap = summer(num_lines)*0.85; 
cmap = cmap(end:-1:1,:); 
int2_ycoord = 0.25; 
src2_ycoord = 0.75; 
xvc_cnt = 1; 
for i = x_loc2plt 
    color2plt = cmap(xvc_cnt,:);
    tgt_psp = reshape(Tgt_PSP(:,i,y_loc2plt), [lenT, 1]); 
    plot(ax_tgt, t, tgt_psp,'LineWidth',1.5,'Color',color2plt);
    
    int2_spkt = Int2_spkt{i,y_loc2plt};
    for spkt_idx = 1:length(int2_spkt)
        plot(ax_spkt, int2_spkt(spkt_idx), int2_ycoord, ... 
            'MarkerEdgeColor', 'w',...
            'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', color2plt, 'color', color2plt);
    end
    
    src2_spkt = Src2_spkt{i,y_loc2plt};
    for spkt_idx = 1:length(src2_spkt)
        plot(ax_spkt, src2_spkt(spkt_idx), src2_ycoord, ...
            'MarkerEdgeColor', 'w', ...
            'Marker', 'd', 'MarkerSize', 8,'MarkerFaceColor', color2plt, 'color', color2plt);
    end
    
    xvc_cnt = xvc_cnt + 1; 
end

text(ax_spkt, src2_spkt(spkt_idx) + 0.5, src2_ycoord + 0.02, 'Src_2', 'Fontsize', 10) ;
text(ax_spkt, int2_spkt(spkt_idx) + 0.5, int2_ycoord + 0.02, 'Int_2', 'Fontsize', 10) ;

colormap(ax_tgt, cmap);
clrbar = colorbar(ax_tgt); 
clrbar_values = x_vec(x_loc2plt); 
tick_pos = linspace(0,1,length(x_loc2plt));
tick_skip = 2; 
clrbar_show = clrbar_values(1:tick_skip:end); 
tick_show = tick_pos(1:tick_skip:end); 
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', tick_show, 'ticklabels', clrbar_show, ...
    'Position', pos_clrbar, 'FontSize', 10);
title(clrbar, '\Deltat_{inp} (ms)','fontsize',11);

xlim(ax_tgt, [70, 80]);
ylim(ax_tgt, [-5,5]); 
tscale_start = 77; 
pspscale_start = 1.5; 
plot(ax_tgt, tscale_start+[0,time_scale], pspscale_start*[1,1], ...
    '-k', 'LineWidth', 0.75); 
plot(ax_tgt, (tscale_start+time_scale)*[1,1], pspscale_start+[0,psp_scale], ...
    '-k', 'LineWidth', 0.75); 
ttl_tgt = title(ax_tgt, 'PSP of Tgt', 'fontweight', 'normal', 'fontsize', title_size); 
set(ax_tgt, 'xcolor', 'none', 'ycolor', 'none'); 
set(ax_tgt, 'visible', 'off'); 
set(ttl_tgt, 'visible', 'on');

xlim(ax_spkt, [70, 80]); 
ylim(ax_spkt, [-0.1,2]);
set(ax_spkt, 'visible', 'off'); 

%% Fig 2C-E
splt_tag = {'C', 'D', 'E'};  
field2plt = {'peak', 'iw', 'auc'}; 

y_val2plt = 6:10 ; 
y_loc2plt = arrayfun(@(rr) find(abs(y_vec-rr)<eps), y_val2plt, 'uniformoutput', true); 

cmap = bone(length(y_loc2plt)+1);
cmap(end,:) = [];
    
for i = 1:length(splt_tag) 
    tag_i = splt_tag{i};
    
    create_ann(['ann_' tag_i], tag_i); 
    
    ax_dat = create_fig(['fig_' tag_i]);
    hold(ax_dat, 'on'); 
    field_i = psp_anly.(field2plt{i}); 
    data2plt_i = field_i.dat; 
    title_i = field_i.name; 
    unit_i = field_i.unit; 
    
    max_yval = 0;  
    min_yval = 100;  
    for j = 1:length(y_loc2plt)
        data2plt_ij = data2plt_i(:,y_loc2plt(j)); 
        plot(ax_dat, x_vec, data2plt_ij, '-' , ...
            'LineWidth', 2 ,'Color',cmap(j,:));
        max_yval = max([max_yval, max(data2plt_ij(:))]); 
        min_yval = min([min_yval, min(data2plt_ij(:))]); 
    end
    title(ax_dat, title_i, 'fontweight', 'normal', 'fontsize', title_size); 
    xlim(ax_dat, [-0.3, 5.1]); 
    ylim(ax_dat, [min_yval, max_yval] + [-1,+1]*(max_yval - min_yval) * 0.1); 
    xlabel(ax_dat, '\Deltat_{inp} (ms)', 'fontsize', title_size-1); 
    ylabel(ax_dat, unit_i); 
    
    if strcmp(tag_i,'E')
        colormap(ax_dat, cmap);
        clrbar = colorbar(ax_dat);
        pos_clrbar = layout_map('colorbar_E').normz_pos;
        
        clrbar_values = y_vec(y_loc2plt);
        tick_pos = linspace(0,1,length(y_loc2plt));
        tick_skip = 1;
        clrbar_show = clrbar_values(1:tick_skip:end);
        tick_show = tick_pos(1:tick_skip:end);
        set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
            'ticks', tick_show, 'ticklabels', clrbar_show, ...
            'Position', pos_clrbar, 'YDir', 'reverse');
        title(clrbar, 'G_{GABA\rightarrowTgt}(nS)');
    end
end
