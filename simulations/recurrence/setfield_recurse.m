function s = setfield_recurse(s, f)
f = strsplit(f, '.'); 
s = setfield(s, f{:});  
end