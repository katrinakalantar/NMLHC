function [fv dfv] = rlrobj(w, gamma, x, y, lambda, options)

% compute regularisation term 
[regV regDV] = feval(options.regFunc, w, options.sn);

t     = x * w;

s0    = (gamma(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (gamma(2,1) ./ (1+exp(-t)));
s1    = (gamma(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (gamma(2,2) ./ (1+exp(-t)));

s0(s0==0) = eps; 
s1(s1==0) = eps;

fv = -sum((repmat(y,1,size(w,2)) .* log(s1)) + (repmat(1-y,1,size(w,2)).* log(s0)),1)...
        + sum(lambda .* regV);

dfv = w(:,1)';
