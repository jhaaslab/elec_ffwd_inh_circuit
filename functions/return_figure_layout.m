function [layout_map, dimensions] = return_figure_layout(file_name)
svg_file = xml2struct(file_name);
file_attr = [svg_file.Attributes];
attr_map = containers.Map({file_attr.Name},{file_attr.Value});
width_str = attr_map('width');
height_str = attr_map('height');
dim_unit = width_str(regexpi(width_str, '[a-z]'));
max_width = str2double(width_str(regexpi(width_str, '\d')));
max_height = str2double(height_str(regexpi(height_str, '\d')));

nonempty_pos = @(c_l) cellfun(@(x) ~isempty(x), c_l, 'UniformOutput', true);
return_nonempty_obj = @(raw_obj) raw_obj(nonempty_pos({raw_obj.Attributes}));
return_named_obj = @(parent, name) return_nonempty_obj([parent(strcmp({parent.Name}, name)).Children]);
return_atrr_val = @(obj, attr) obj.Attributes(strcmp({obj.Attributes.Name},attr)).Value;
file_children = return_nonempty_obj([svg_file.Children]);
graphic_obj = return_named_obj(file_children, 'g');
for i = 1:length(graphic_obj)
    if strcmp(return_atrr_val(graphic_obj(i),'label'), 'layout')
        layout_obj = return_nonempty_obj([graphic_obj(i).Children]);
        break
    end
end
layout_map = containers.Map();
layout_params = {'x', 'y', 'width', 'height'};
normz_fun = {@(x) x/max_width, @(x) 1-(x/max_height), @(x) x/max_width, @(x) x/max_height};

for i = 1:length(layout_obj)
    ly_i = layout_obj(i);
    lbl = return_atrr_val(ly_i, 'label');
    tmp_struct = struct();
    for j = 1:length(layout_params)
        normz_j = normz_fun{j};
        val_j = normz_j(str2double(return_atrr_val(ly_i, layout_params{j}))); 
        tmp_struct.(layout_params{j}) = val_j;
    end
    tmp_struct.y = tmp_struct.y - tmp_struct.height; 
    tmp_struct.normz_pos = [tmp_struct.x, tmp_struct.y, tmp_struct.width, tmp_struct.height];
    layout_map(lbl) = tmp_struct;
end
dimensions = struct('width', max_width, 'height', max_height, 'unit', dim_unit); 
end