function s = setfield_recurse(s, f, v)
f = strsplit(f, '.'); 
s = setfield(s, f{:}, v);  
end