function [v dv] = l2(w, void)

v  = sum(w.^2); 
dv = 2*w;