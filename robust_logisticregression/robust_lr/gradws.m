% function definition for 'minFunc' optimiser
% fv  = function value
% dfv = gradient of the fucntion w.r.t W

function [fv dv] = gradws(w, gamma, x, y, mask, lambda, options)

CLS    = size(gamma, 2);
DIM    = size(x, 2);

% compute regularisation term 
[regV regDV] = feval(options.regFunc, w, options.sn);

reg    = sum(sum(lambda.*regV));

w      = reshape(w, DIM, CLS);
lambda = reshape(lambda, DIM, CLS);
regDV  = reshape(regDV, DIM, CLS);

% calculate function value
sk     = posteriorYHat(x, w, gamma);    
fv     = -sum(log(sk(mask)), 1) + reg;

% faster way
tmpw   = zeros(DIM, CLS); 
logit  = softmax(x * w);

for c = 1:CLS 
    
    % calculate logit w.r.t to w_c and repmat to match sk
    logit_c = repmat(logit(:,c), 1, CLS);
    
    % calculate the big chunk in parentheses
    gc      = repmat(gamma(c,:), CLS, 1);
    gc      = gc - gamma;    
    chunk   = logit * gc;
    
    % calculate multiplier to X
    multi   = (logit_c .* chunk) ./ sk ;
    
    % multiply with X
    llhw    = zeros(1, DIM); 
    
    for k = 1:CLS
        %meaning res   = repmat(multi(:,k), 1, DIM) .* x;
        res   = bsxfun(@times,multi(:,k),x);
        llhw  = llhw - sum(res(y==k,:), 1);
    end
    
    % save the derivative of w_c
    tmpw(:,c) = llhw';
end

% compute the derivative of the objective
tmpw = tmpw + (lambda .* regDV);
dv   = cat(2, tmpw(:));

 