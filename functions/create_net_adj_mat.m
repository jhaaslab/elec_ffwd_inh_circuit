clc; clear;
Nv = [5,10,20,50];
fid = fopen('prm_files\extended_net_con.yaml', 'w');
for i = 1:length(Nv)
    net_str = create_net_str(Nv(i),i-1); 
    fprintf(fid, net_str);
end
fclose(fid);

function net_str = create_net_str(N, net_id)
net_str = sprintf('Net_%d:\n', net_id);
node_str = sprintf(' nodes:\n');
list_nodes = {};
for i = 1 : N
    node_str = sprintf('%s  - [''%s_%02d'', ''%s'']\n', ...
        node_str, 'Src', i, 'RS');
    node_str = sprintf('%s  - [''%s_%02d'', ''%s'']\n', ...
        node_str, 'Int', i, 'FS');
    node_str = sprintf('%s  - [''%s_%02d'', ''%s'']\n', ...
        node_str, 'Tgt', i, 'RS');
    list_nodes = [list_nodes, ...
        {'Src', 'Int', 'Tgt';...
        i, i, i}]; %#ok<AGROW>
end
edge_str = sprintf(' edges:\n');

for i = 1 : (N*3)
    node_i = list_nodes{1,i};
    id_i = list_nodes{2,i};
    edge_i = sprintf('  - [');
    for j = 1 : (N*3)
        node_j = list_nodes{1,j};
        id_j = list_nodes{2,j};
        if id_i == id_j
            switch upper(node_i)
                case 'SRC'
                    if strcmpi(node_j,'INT') || strcmpi(node_j,'TGT')
                        edge_i = sprintf('%sAMPA',edge_i);
                    else
                        edge_i = sprintf('%s0',edge_i);
                    end
                case 'INT'
                    if strcmpi(node_j,'TGT')
                        edge_i = sprintf('%sGABA',edge_i);
                    else
                        edge_i = sprintf('%s0',edge_i);
                    end
                otherwise
                    edge_i = sprintf('%s0',edge_i);
            end
        else
            if strcmpi(node_i,'INT') && strcmpi(node_j,'INT')
                edge_i = sprintf('%sELEC-GABA',edge_i);
            else
                edge_i = sprintf('%s0',edge_i);
            end
        end
        if j < N*3
            edge_i = sprintf('%s,',edge_i);
        else
            edge_i = sprintf('%s]',edge_i);
        end
    end
    edge_str = sprintf('%s%s\n',edge_str,edge_i);
end
net_str = sprintf('%s%s%s\n',net_str,node_str,edge_str);
end
