function [setup_struct, simvar_struct, field_names, var_alias] = find_sim_var(setup_struct, field_list)
import setup_sweep.* 

ind_with_vec = ~cellfun(@isempty, regexp(field_list, '.vec$', 'match'), 'uni', 1);
vec_fields = field_list(ind_with_vec);
var_vecs = cellfun(@(x) getfield_recurse(setup_struct, x), vec_fields, 'uni', 0);

field_names = cellfun(@(x) x(1:end-4), vec_fields, 'uni', 0);

lbl_fields = cellfun(@(x) [x, '.label'], field_names, 'uni', 0);
has_lbl = find(cellfun(@(x) any(strcmp(field_list, x)), lbl_fields, 'uni', 1));

var_labels = field_names;
var_labels(has_lbl) = arrayfun(@(x) getfield_recurse(setup_struct, lbl_fields{x}), has_lbl, 'uni', 0);

var_alias = regexprep(field_names, '\.', '_');

nonsim_fields = field_list(~ind_with_vec); 
for i = 1:length(nonsim_fields) 
    val_i = getfield_recurse(setup_struct, nonsim_fields{i}); 
    if iscell(val_i)
        if isnumeric(val_i{1}) 
            setup_struct = setfield_recurse(setup_struct, ...
                nonsim_fields{i}, cell2mat(val_i)); 
        end
    end
end


simvar_struct = struct(); 
for i = 1:length(field_names)
    vec_var_i = var_vecs{i}; 
    if iscell(vec_var_i{1}) 
        if isnumeric(vec_var_i{1}{1})
            vec_var_i = cellfun(@cell2mat, vec_var_i, 'uni', 0);
        end
    end
    setup_struct = setfield_recurse(setup_struct, field_names{i}, []); 
    simvar_struct.(var_alias{i}) = struct('vec', {vec_var_i}, 'label', var_labels{i}); 
end

end
