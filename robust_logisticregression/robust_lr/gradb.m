% function definition for Rasmussen's minimise()
% fv  = function value
% dfv = gradient of the function w.r.t W

function [fv dv] = gradb(sqb, w, K, y, g, reg_w, reg_b)

% compute regularisation term 
reg   = 0.5 * reg_w * sum(w(2:end).^2) + sum(reg_b .* (sqb.^2)); 

x     = compositeKernel(sqb.^2, K);
t     = x * w;

s0    = (g(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,1) ./ (1+exp(-t)));
s1    = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));

s0(s0==0) = eps; 
s1(s1==0) = eps;

fv    = -sum((y .* log(s1)) + ((1 - y).* log(s0)),1) + reg;

% Eq.34 without product with x_n
tmp1  = (g(2,2) - g(1,2)) * (y ./ s1);      
tmp2  = (g(2,1) - g(1,1)) * ((1 - y) ./ s0); 

gAux1 = (tmp1 + tmp2) ./ (1+exp(-t)) .* (1 ./ ((1 ./ exp(-t))+1));

len   = length(sqb);
dv    = zeros(len,1);
reg   = reg_b .* 2 .* sqb; 
for i=1:len    
    wt    = w * 2 * sqb(i);    
    dv(i) = -sum(gAux1 .* (K{i} * wt), 1) + reg(i);    
end

