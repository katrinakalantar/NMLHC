function [v dv] = lasso(w, sn)

v  = sqrt(w.^2 + sn);
dv = w ./ (sqrt(w.^2 + sn)); 