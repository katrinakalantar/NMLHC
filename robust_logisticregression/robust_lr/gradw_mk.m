% function definition for minFunc()
% fv  = function value
% dfv = gradient of the function w.r.t W

function [fv dfv] = gradw_mk(w, x, y, g, b, reg_w, reg_b)

% compute regularisation term 
reg  = 0.5 * reg_w * sum(w(2:end).^2) +  sum(reg_b .* b); 

t    = x * w;

s0   = (g(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,1) ./ (1+exp(-t)));
s1   = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));

s0(s0==0) = eps; 
s1(s1==0) = eps;

fv   = -sum((y .* log(s1)) + ((1 - y).* log(s0)),1) + reg;

% compute derivative of regularisation term 
reg  = [0; reg_w * w(2:end)];

% Eq.34 without product with x_n
tmp1  = (g(2,2) - g(1,2)) * (y ./ s1);      
tmp2  = (g(2,1) - g(1,1)) * ((1 - y) ./ s0); 

gAux1 = (tmp1 + tmp2) ./ (1+exp(-t)) .* (1 ./ ((1 ./ exp(-t))+1));

% completed Eq.34 with dot product with x
%gAux2 = repmat(gAux1, 1, size(x,2)) .* x;  % slower

gAux2 = bsxfun(@times, gAux1, x); % faster

dfv   = -sum(gAux2,1)' + reg;
