function lst = recursive_fieldnames(s)
[lst, cnt, ~] = recurse_struct(s); 
lst = lst(1:cnt-1); 
end

function [lst, cnt, imax] = recurse_struct(s, nme, lst, cnt, imax)
max_size = 1e3; 
if nargin == 1
    lst = cell(max_size, 1); 
    cnt = 1; 
    imax = max_size; 
    nme = '';
end

if isstruct(s)
    
    if ~isempty(nme) 
        nme = [nme '.'];
    end
    
    fns = fieldnames(s);
    for i = 1:length(fns)
        nme_i = [nme fns{i}];
        strc_i = s.(fns{i});
        [lst, cnt, imax] = recurse_struct(strc_i, nme_i, lst, cnt, imax);
    end
    
else
    
    if cnt > imax
        imax = imax + max_size; 
        lst{imax} = ''; 
    end
    lst{cnt} = nme;
    cnt = cnt + 1; 
    
end

end