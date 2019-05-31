function tmp_rect = template_rectanglepulse(T, fact_dur)
lenT = length(T); 

tmp_rect = zeros(1, lenT); 
dur_rect = ceil(fact_dur*lenT);
tmp_rect(1:dur_rect) = 1;
end
