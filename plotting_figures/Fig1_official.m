clc; clear; close all; 
addpath(genpath('../figures'), genpath('../functions')); 
%% Load data
load('../data/net0_sim/Net0_AMPAvsGABA.mat')
x_name = data2save.vars(1,:);
x_vec = data2save.(x_name);
y_name = data2save.vars(2,:);
y_vec = data2save.(y_name);

[xz, yz] = meshgrid(x_vec,y_vec);
avail_attr = {'peak_epsp', 'dt_epsp', 'auc_epsp'};
attr_units = {'mV', 'ms', 'ms x mV'};
valid_plot_name = @(nme) [nme(1:regexp(nme, '_')) '{' nme(regexp(nme, '_')+1:end) '}'];
labels_ = struct();
labels_.x = 'G_{AMPA\rightarrowTgt} (nS)';
labels_.y = 'G_{GABA\rightarrowTgt} (nS)';
labels_.title = {'Peak of PSP', 'Integration window', 'AUC of PSP'};
labels_.colorbar = {'mV', 'ms', 'ms\cdotmV'};

rnd_dat = randi(data2save.num_dat);
dat_ = data2save.(sprintf('dat_%06d',rnd_dat));
t = data2save.dat_000001.dat.t;
Int_spkt = zeros(size(xz));
Src_spkt = zeros(size(xz));
Tgt_PSP = zeros([length(t), size(xz)]);


for i = 1:data2save.num_dat
    nme_ = sprintf('dat_%06d',i); 
    dat_ = data2save.(nme_);
    data2save.(nme_).iw = analyze_TgtPSP(dat_.dat.t, dat_.dat.Tgt);
end

for i = 1:data2save.num_dat
    dat_ = data2save.(sprintf('dat_%06d',i));
    var_ = dat_.pair;
    x_loc = find(abs(x_vec-var_(1))<eps);
    y_loc = find(abs(y_vec-var_(2))<eps);
    [~,spkt] = findpeaks(dat_.dat.Src, t, 'MinPeakProminence', 50);
    Src_spkt(x_loc,y_loc) = spkt;
    [~,spkt] = findpeaks(dat_.dat.Int, t, 'MinPeakProminence', 50);
    Int_spkt(x_loc,y_loc) = spkt;
    Tgt_PSP(:,x_loc,y_loc) = dat_.dat.Tgt;
end

%% Load layout
file_name  = 'Fig1_layout_official.svg';
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
ann_style = {'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontSize', 15, 'FontWeight', 'bold'}; 

create_ann = @(tag_name, string_) annotation('textbox', 'Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos, ...
    'String', string_, ann_style{:}); 
create_fig = @(tag_name) axes('Units', 'normalized', ...
    'Position', layout_map(tag_name).normz_pos);

title_size = 10; 
%% Fig 1A1
create_ann('ann_A', 'A'); 
for i=1:4
    ann_name = ['ann_A' num2str(i)];    
    create_ann(ann_name, num2str(i)); 
end

%% Fig 1B
create_ann('ann_B', 'B');

create_fig('fig_B');

hold on;
plot([0,100], [0,0], ':', 'color', 0.1*[1,1,1], 'linewidth',0.75);
tgt_ = Tgt_PSP(:,3,8);
tgt_ = tgt_ - tgt_(1);
plot(t, tgt_, 'LineWidth', 1.5, 'Color','k');
loc_iw = find(tgt_ > 0);
for i=round(linspace(2,length(loc_iw)-1,10))
    plot([t(loc_iw(i)), t(loc_iw(i))+0.1], [0, tgt_(loc_iw(i))], '-k');
end
plot([t(loc_iw(2)), t(loc_iw(end-1))], [-.5,-.5], '.-k', 'MarkerSize', 9);
plot((t(loc_iw(end))+1)*[1,1], [0.3, max(tgt_)], '.-k', 'MarkerSize', 9);
notation_size = 8; 
text(t(loc_iw(3))-1.1, -1.6, {'Integration','window'}, 'Fontsize', notation_size);
text(t(loc_iw(end))+1.5, 1.10, 'Peak_{PSP}', 'Fontsize', notation_size);
text(27.5, 1.5, 'AUC_{PSP}', 'Fontsize', notation_size);
ylim([-4,3]); xlim([28,45]);
set(gca, 'visible', 'off');
title('PSP of Tgt', 'fontweight', 'normal', 'fontsize', title_size, 'visible', 'on'); 
%% Fig 1C
g_ampa2look = 3;
g_gaba2look = 6;

create_ann('ann_C', 'C');

create_fig('fig_C');

hold on;
cmap = bone(length(y_vec)+1); 
cmap(end,:) = ''; 
plot([0,100], [0,0], ':', 'color', 0.6*[1,1,1], 'linewidth',0.75);
for i = 1:length(y_vec)
    x_loc = find(abs(x_vec-g_ampa2look)<eps);
    y_loc = find(abs(y_vec-y_vec(i))<eps);
    tgt_i = Tgt_PSP(:,x_loc,y_loc);
    plot(t, tgt_i - tgt_i(1), 'LineWidth', 1.5, 'Color', cmap(i,:));
end
src_i = Src_spkt(x_loc, y_loc);
int_i = Int_spkt(x_loc, y_loc);
src_plt = plot([src_i, src_i], [-5, 5], '-k', 'DisplayName', 'Src', 'LineWidth', 0.75); %#ok<NASGU>
text(src_i - 6, -4.5, 'Src', 'fontsize', 10); 
int_plt = plot([int_i, int_i], [-5, 5], '--k', 'DisplayName', 'Int', 'LineWidth', 0.75); %#ok<NASGU>
text(int_i + 1.5, -4.5, 'Int', 'fontsize', 10); 
plot([50, 55], [-4, -4], '-k', 'LineWidth', 0.75);
plot([55, 55], [-4, -3], '-k', 'LineWidth', 0.75);
ylim([-5,5]); xlim([28,60]);
clrbar = colorbar(gca); colormap(gca, cmap); caxis([min(y_vec), max(y_vec)]);
set(clrbar, 'fontsize', 8, 'Box', 'off', 'LineWidth', 0.001);
ylabel(clrbar, labels_.y, 'fontsize', 10);
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'visible', 'off')
title(gca, 'PSP of Tgt', ...
    'fontweight', 'normal', 'fontsize', title_size, 'visible', 'on');
%% Fig 1D
create_ann('ann_D', 'D');

create_fig('fig_D');

hold on;
cmap = bone(length(y_vec)+1); 
cmap(end,:) = ''; 

plot([0,100], [0,0], ':', 'color', 0.6*[1,1,1], 'linewidth',0.75);
for i = 1:length(x_vec)
    x_loc = find(abs(x_vec-x_vec(i))<eps);
    y_loc = find(abs(y_vec-g_gaba2look)<eps);
    tgt_i = Tgt_PSP(:,x_loc,y_loc);
    plot(t, tgt_i - tgt_i(1), 'LineWidth', 1.5, 'Color', cmap(i,:));
    
end
src_i = Src_spkt(x_loc, y_loc);
int_i = Int_spkt(x_loc, y_loc);
src_plt = plot([src_i, src_i], [-3, 8], '-k', 'DisplayName', 'Src', 'LineWidth', 0.75);
int_plt = plot([int_i, int_i], [-3, 8], '--k', 'DisplayName', 'Int', 'LineWidth', 0.75);
plot([50, 55], [0.5,0.5], '-k', 'LineWidth', 0.75);
plot([55, 55], [0.5,1.5], '-k', 'LineWidth', 0.75);
ylim([-3,7]); xlim([28,60]);
clrbar = colorbar(gca); colormap(gca,cmap); caxis([min(x_vec), max(x_vec)]);
set(clrbar, 'fontsize', 8, 'Box', 'off', 'LineWidth', 0.001);
ylabel(clrbar, labels_.x, 'fontsize', 10);
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'visible', 'off')
title(gca, 'PSP of Tgt', ...
    'fontweight', 'normal', 'fontsize', title_size, 'visible', 'on');
%% Fig 1E-G
data2plt = zeros(size(xz));
splt_names = {'E', 'F', 'G'};
for j = 1:length(avail_attr)
    attr2plt = avail_attr{j};
    unit2plt = attr_units{j};
    
    for i = 1:data2save.num_dat
        dat_ = data2save.(sprintf('dat_%06d',i));
        var_ = dat_.pair;
        x_loc = find(abs(x_vec-var_(1))<eps);
        y_loc = find(abs(y_vec-var_(2))<eps);
        data2plt(x_loc,y_loc) = dat_.iw.(attr2plt);
    end
    create_ann(['ann_' splt_names{j}], splt_names{j});
    create_fig(['fig_', splt_names{j}]);
    image(x_vec,y_vec,data2plt', 'CDataMapping','scaled');
    clrbr = colorbar(gca);
    title(clrbr,labels_.colorbar{j},'fontsize', 7);
    set(clrbr, 'Box', 'off')
    colormap(gca,'hot')
    cmap = colormap(gca);
    cmap(end,:) = ''; 
    colormap(gca,cmap); 
    xlabel(labels_.x, 'fontsize', 10);
    ylabel(labels_.y, 'fontsize', 10);
    title(labels_.title{j}, 'fontweight', 'normal', 'fontsize', title_size);
    set(gca, 'ydir', 'normal')
    daspect([1,1,1]);
end
