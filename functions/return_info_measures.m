function info_res =  return_info_measures(data_folder, filename_prefix, ...
    saved_name, cell_pref, bw_opt, smw_opt)
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
value_var = cellfun(@(x) double(vars_.(x).range), var_names,'UniformOutput',false);
len_var = cellfun(@(x) length(x), value_var, 'UniformOutput', true);

cell_names = cellstr(data2save.cell_names)';

num_pref = length(cell_pref);
idx_pref  = struct();
for i = 1:num_pref
    idx_pref.(cell_pref{i}) = find(contains(cell_names, cell_pref{i}));
end
raw_data = cell([len_var',num_pref]);
raw_co_spike = cell(len_var'); 

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
    
    num_unit = length(idx_pref.(cell_pref{1}));
    co_occ = nan(num_pref,num_unit);
    for idx_i = 1:num_unit
        for j = 1:num_pref
            idx2find = find(strcmp(cell_names,sprintf('%s_%02d',cell_pref{j},idx_i)));
            spkt_j = spkts(index == idx2find);
            if ~isempty(spkt_j) 
                co_occ(j,idx_i) = spkt_j;
            end
            if length(spkt_j) > 1
                error('The combo [%s], at %s_%02d has more than 1 spike', ...
                    sprintf(' %d ',[loc_{:}]),cell_pref{j},idx_i);
            end
        end       
    end
    raw_co_spike{loc_{:}} = co_occ;
        
end

sizeRD = size(raw_data);
stat_dim = num2cell(sizeRD(2:end)); 
stat_dist(stat_dim{:}) = struct(...
    'mi_src_out', [], 'mi_over_hsrc', [], 'mi_over_hout', [], 'entropy_lat_out2src', []); 

grid_arrays = arrayfun(@(x) {1:x}, sizeRD(2:end-1)); 
grid_cells = cell(1,length(grid_arrays));
[grid_cells{:}] = ndgrid(grid_arrays{:}); 
tot_combos = numel(grid_cells{1}(:)); 
grid_combos = zeros(tot_combos, length(grid_arrays)); 
for i = 1:length(grid_arrays)
    grid_combos(:,i) = grid_cells{i}(:);
end
Src_loc = find(strcmp(cell_pref, 'Src')); 
t_init_std = value_var{2};


smw_calc_mi = smw_opt;
sm_1D = hanning(smw_calc_mi);
sm_1D = sm_1D/sum(sm_1D);
sm_2D = sm_1D*sm_1D';
sm_2D = sm_2D/sum(sm_2D(:));
inf_funcs = information_analysis_functions;

if ischar(bw_opt)
    switch lower(bw_opt)
        case 'fixed'
            f_binw_calc_mi = @(x) (0.01);
        case 'variable'
            f_binw_calc_mi = @(x) (x/100);
        otherwise
            error(['return_info_measures::''bw_opt'' can only be either ''fixed'' or ''variable''' ...
                'if it''s a string, not %s'], bw_opt);
    end
else
    if isa(bw_opt, 'function_handle')
        f_binw_calc_mi = bw_opt; 
    else
        error('return_info_measures::''bw_opt'' can only be either be a string or a function handle');
    end
end

for i = 1:tot_combos
    combo_i = num2cell(grid_combos(i,:)); 
    t_init_std_i = t_init_std(combo_i{1}); 
    Src_dat = [raw_data{:,combo_i{:},Src_loc}];     
    mean_src = mean(Src_dat); 
    
    raw_cospk_i = [raw_co_spike{:,combo_i{:}}];
    src_col = raw_cospk_i(Src_loc,:) - mean_src; 
    binw_calc_mi = f_binw_calc_mi(t_init_std_i); 

    for j = 1:sizeRD(end) 
        ij_col = raw_cospk_i(j,:)-mean_src; 
        info_res = inf_funcs.return_smoothed_hist_2D_and_MI(src_col,1,ij_col,1,sm_1D,sm_2D);
        stat_dist(combo_i{:},j).mi_src_out = info_res.MI;
        stat_dist(combo_i{:},j).mi_over_hsrc = info_res.MI/info_res.dat1.entropy; 
        stat_dist(combo_i{:},j).mi_over_hout = info_res.MI/info_res.dat2.entropy; 
        
        lat_out2src = ij_col - src_col; 
        info_lat = inf_funcs.return_smoothed_hist_1D(lat_out2src,binw_calc_mi,sm_1D);
        stat_dist(combo_i{:},j).entropy_lat_out2src = info_lat.entropy;
    end

end

info_res = struct(); 
info_res.info_measures = stat_dist; 
info_res.var_names = var_names; 
info_res.value_var = value_var; 

if ischar(saved_name) 
    save([saved_name '.mat'], 'info_res');
end
end