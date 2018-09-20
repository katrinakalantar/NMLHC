% function definition for 'minFunc' optimiser
% fv  = function value
% dfv = gradient of the function w.r.t W

function [fv dfv] = gradw(w, gamma, x, y, lambda, options)

% compute regularisation term 
[regV regDV] = feval(options.regFunc, w, options.sn);

t     = x * w;

s0    = (gamma(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (gamma(2,1) ./ (1+exp(-t)));
s1    = (gamma(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (gamma(2,2) ./ (1+exp(-t)));

s0(s0==0) = eps; 
s1(s1==0) = eps;

% set the distribution of data 
D     = options.dist;   
fv    = -sum((D .* y .* log(s1)) + (D .* (1 - y).* log(s0)),1) + sum(lambda .* regV);
% compute corresponding derivative
tmp1  = (gamma(2,2) - gamma(1,2)) * (D .* y ./ s1);
tmp2  = (gamma(2,1) - gamma(1,1)) * (D .* (1 - y) ./ s0);

       
gAux1 = (tmp1 + tmp2) ./ (1+exp(-t)) .* (1 ./ ((1 ./ exp(-t))+1));

% completed Eq.34 with dot product with x_n
gAux2 = repmat(gAux1, 1, size(x,2)) .* x;
dfv   = -sum(gAux2,1)' + (lambda .* regDV);

