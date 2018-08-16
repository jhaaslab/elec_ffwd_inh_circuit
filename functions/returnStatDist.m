function [stat_dist, x_axis_actual, y_axis_actual] = returnStatDist(data_folder, filename_prefix, cell_pref)

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
% unit_var = cellfun(@(x) vars_.(x).unit, var_names,'UniformOutput',false);
value_var = cellfun(@(x) double(vars_.(x).range), var_names,'UniformOutput',false);
len_var = cellfun(@(x) length(x), value_var, 'UniformOutput', true);

cell_names = cellstr(data2save.cell_names)';
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
    Src_dat = [raw_data{:,combo_i{:},Src_loc}]; %#ok<FNDSB>
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
        stat_dist(combo_i{:},j).lat10per = t_10percent_srcpeak/std_src;
    end
    tinp_dist{combo_i{:}} = [t_inp{:,combo_i{:}}];
end

%% Value vectors
t_init_std = value_var{2};
g_elec_vec = value_var{3};

f_vec2show = @(x) [x(1), x(ceil(length(x)/2)), x(end)];
x_axis_actual = struct(...
    'name', '\sigma_{inp} (ms)', ...
    'vec', t_init_std, ...
    'show', f_vec2show(t_init_std) );
y_axis_actual = struct(...
    'name', '\SigmaG_{elec} (nS)', ...
    'vec', g_elec_vec, ...
    'show', f_vec2show(g_elec_vec) );

end