function [raw_co_spike, num_unit_trial, sigma_inp, g_elec] = return_co_spike(data_folder, filename_prefix) 
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
cell_pref = {'Src', 'Int', 'Tgt'};
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
    % no spike, for simplicity of visualization
    % this number is arbitrary, as long as it is outside of the range of
    % possible latencies, when calculating entropy, as it can be
    % considered as anothre state. For actual calculation of entropy and
    % mutal information (refer to "information_analysis_functions"), the 
    % value to replace NaN (no spike) is 2 standard deviations from the
    % maximum possible values of the distribution of interest. Again, this
    % is arbitrary, as long as the replacement value is outside of the
    % possible values of the distribution. 
    co_occ(isnan(co_occ)) = 1000; 
    raw_co_spike{loc_{:}} = co_occ;
end

num_trials = len_var(order_var(strcmp(var_names, 'n_trial')));
num_unit = length(idx_pref.(cell_pref{1}));
num_unit_trial = num_trials * num_unit;

sigma_inp = [value_var{order_var(strcmp(var_names, 't_init_std'))}];
g_elec = [value_var{order_var(strcmp(var_names, 'g_elec_int_to_int'))}];


