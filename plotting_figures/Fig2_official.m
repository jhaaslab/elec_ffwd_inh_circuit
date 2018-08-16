clc; clear; close all; 
addpath(genpath('../figures'), genpath('../functions')); 
%% Load data 
load('../data/net1_sim/Net1_dtInpvsholdInt.mat')

for i = 1:data2save.num_dat
    nme_ = sprintf('dat_%06d',i); 
    dat_ = data2save.(nme_);
    data2save.(nme_).iw = analyze_TgtPSP(dat_.dat.t, dat_.dat.Tgt);
end

var_names = fieldnames(data2save.vars);
x_name = var_names{1};
x_vec = data2save.vars.(x_name);
y_name = var_names{2};
y_vec = data2save.vars.(y_name);

valid_plot_name = @(nme) [nme(1:regexp(nme, '_')) '{' nme(regexp(nme, '_')+1:end) '}'];

t = data2save.dat_000001.dat.t;
lenT = length(t);

Int_spkt = cell(length(x_vec),length(y_vec));
Src1_spkt = []; 
Src2_spkt = cell(length(x_vec),length(y_vec));
Tgt_PSP = zeros(lenT,length(x_vec),length(y_vec));
iw_psp = zeros(length(x_vec),length(y_vec));
peak_psp = zeros(length(x_vec),length(y_vec));
auc_psp = zeros(length(x_vec),length(y_vec));
for i = 1:data2save.num_dat
    dat_ = data2save.(sprintf('dat_%06d',i));
    var_ = dat_.var;
    x_loc = find(abs(x_vec-var_(1))<eps);
    y_loc = find(abs(y_vec-var_(2))<eps);
    if isempty(x_loc) && isempty(y_loc)
        error('Could not find ''index'' of %s', sprintf('dat_%06d',i));
    end
    tgt_vm = dat_.dat.Tgt; 
    Tgt_PSP(:,x_loc,y_loc) = tgt_vm - tgt_vm(1); 
    
    [~,loc_spk] = findpeaks(dat_.dat.Int, t, 'MinPeakProminence', 50); 
    Int_spkt{x_loc,y_loc} = loc_spk; 
    
    [~,loc_spk] = findpeaks(dat_.dat.Src_2, t, 'MinPeakProminence', 50); 
    Src2_spkt{x_loc,y_loc} = loc_spk;
    
    if isempty(Src1_spkt)
        [~,loc_spk] = findpeaks(dat_.dat.Src_1, t, 'MinPeakProminence', 50);
        Src1_spkt = loc_spk;
    end
    
    iw_psp(x_loc,y_loc) = dat_.iw.dt_epsp;
    peak_psp(x_loc,y_loc) = dat_.iw.peak_epsp;
    auc_psp(x_loc,y_loc) = dat_.iw.auc_epsp;    
end
psp_anly = struct(  'iw', struct('dat', iw_psp, 'name', 'Integration window', 'unit', 'ms'), ...
                    'peak', struct('dat', peak_psp, 'name', 'Peak of PSP', 'unit', 'mV'), ...
                    'auc', struct('dat', auc_psp, 'name', 'AUC of PSP', 'unit', 'ms\cdotmV')); 
clear iw_psp peak_psp auc_psp;

%% Load layout
file_name  = 'Fig2_layout_vtest.svg';
[layout_map, dimensions] = return_figure_layout(file_name);
width = dimensions.width;
height = dimensions.height;
unit = dimensions.unit;
conv_factor = double(unitConversionFactor(str2symunit(unit), str2symunit('cm')));
layout_keys = layout_map.keys();
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height]*conv_factor);

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


ax_spkt = create_fig('fig_B_2');
hold(ax_spkt, 'on'); 
plot(ax_spkt, Src1_spkt*[1,1], [-2.5,3], '--k', 'LineWidth', 0.9); 
text(ax_spkt, Src1_spkt-0.35, -0.45, 'Src_1', 'Fontsize', 10) ;
loc50pA = find(abs(y_vec - 50) < eps); 

ax_tgt = create_fig('fig_B_1'); 
hold(ax_tgt, 'on'); 
plot(ax_tgt, [0, 200], [0, 0], ':', 'color', 0.6*[1,1,1], 'linewidth',0.75);

cond_x_vec = [0,4];
x_vec_cond = x_vec(x_vec >= cond_x_vec(1) & x_vec <= cond_x_vec(2));
x_vec_cond = x_vec_cond(1:4:end);
num_lines = length(x_vec_cond); 
cmap = jet(num_lines);
xvc_cnt = 1; 
for i = 1:length(x_vec)
    if ~isempty(x_vec_cond(x_vec_cond == x_vec(i)))
        color_2plt = cmap(xvc_cnt,:); 
        tgt_psp = reshape(Tgt_PSP(:,i,loc50pA),[lenT,1]);
        plot(ax_tgt, t, tgt_psp,'LineWidth',1,'Color',color_2plt);
        
        int_spkt = Int_spkt{i,loc50pA}; 
        for spkt_idx = 1:length(int_spkt) 
            plot(ax_spkt, int_spkt(spkt_idx), 0.9*xvc_cnt/num_lines, ...
                'Marker', 'o', 'MarkerSize', 3.25, 'MarkerFaceColor', color_2plt, 'color', color_2plt);
        end
        
        src2_spkt = Src2_spkt{i,loc50pA};
        for spkt_idx = 1:length(src2_spkt) 
            plot(ax_spkt, src2_spkt(spkt_idx), 1.2+0.9*xvc_cnt/num_lines, ...
                'Marker', 'd', 'MarkerSize', 3.5,'MarkerFaceColor', color_2plt, 'color', color_2plt);
        end    
        xvc_cnt = xvc_cnt + 1;
    end
end

text(ax_spkt, src2_spkt(spkt_idx) - 1, 1.5, 'Src_2', 'Fontsize', 10) ;
text(ax_spkt, int_spkt(spkt_idx) + 0.5, 0.5, 'Int', 'Fontsize', 10) ;

colormap(ax_tgt, cmap);
clrbar = colorbar(ax_tgt); 
clrbar_show = x_vec_cond(1:5:end);
pos_clrbar = layout_map('colorbar_B').normz_pos;
set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
    'ticks', linspace(0,1,length(clrbar_show)), 'ticklabels', clrbar_show, ...
    'Position', pos_clrbar);
title(clrbar, '\Deltat_{inp} (ms)','fontsize',11);

xlim(ax_tgt, [70, 80]);
ylim(ax_tgt, [-2.5, 4]); 
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
ylim(ax_spkt, [-0.1,2.35]);
set(ax_spkt, 'visible', 'off'); 
ihold_txt = 'I^{hold}_{Int} = 50 pA'; 
text(ax_spkt, 80-3, 0.15, ihold_txt, 'fontsize', 10);

%% Fig 2C-E
splt_tag = {'C', 'D', 'E'};  
field2plt = {'peak', 'iw', 'auc'}; 

for i = 1:length(splt_tag) 
    tag_i = splt_tag{i};
    
    create_ann(['ann_' tag_i], tag_i); 
    
    ax_dat = create_fig(['fig_' tag_i]);
    hold(ax_dat, 'on'); 
    cmap = bone(length(y_vec)+1);
    cmap(end,:) = []; 
    
    field_i = psp_anly.(field2plt{i}); 
    data2plt_i = field_i.dat; 
    title_i = field_i.name; 
    unit_i = field_i.unit; 
    for j = 1:length(y_vec)
        plot(ax_dat, x_vec, data2plt_i(:,j), '-' ,'LineWidth',2,'Color',cmap(j,:));
    end
    title(ax_dat, title_i, 'fontweight', 'normal', 'fontsize', title_size); 
    xlim(ax_dat, [0, 7]); 
    xlabel(ax_dat, '\Deltat_{inp} (ms)', 'fontsize', title_size-1); 
    ylabel(ax_dat, unit_i); 
    if strcmp(tag_i,'E') 
        colormap(ax_dat, cmap);
        clrbar = colorbar(ax_dat);
        pos_clrbar = layout_map('colorbar_E').normz_pos;
        set(clrbar, 'Box', 'off', 'LineWidth', 0.001,...
            'ticks', linspace(0,1,length(y_vec)), 'ticklabels', y_vec, ...
            'Position', pos_clrbar, 'YDir', 'reverse');
        title(clrbar, 'I^{hold}_{Int} (pA)');
    end
    
end
