% function definition for 'minFunc' optimiser
% fv  = function value
% dfv = gradient of the fucntion w.r.t W

function [fv dv] = gradws_mk(w, g, x, y, mask, b, reg_w, reg_b)

CLS    = size(g, 2);
DIM    = size(x, 2);

w      = reshape(w, DIM, CLS);
reg_w  = reshape(reg_w, 1, CLS);

% compute regularisation term 
reg    = 0.5 * (sum(reg_w .* sum(w(2:end,:))).^2) +  sum(reg_b .* b); 


% calculate function value
sk     = posteriorYHat(x, w, g);    
fv     = -sum(log(sk(mask)), 1) + reg;

% faster way
tmpw   = zeros(DIM, CLS); 
logit  = softmax(x * w);

for c = 1:CLS 
    
    % calculate logit w.r.t to w_c and repmat to match sk
    logit_c = repmat(logit(:,c), 1, CLS);
    
    % calculate the big chunk in parentheses
    gc      = repmat(g(c,:), CLS, 1);
    gc      = gc - g;    
    chunk   = logit * gc;
    
    % calculate multiplier to X
    multi   = (logit_c .* chunk) ./ sk ;
    
    % multiply with X
    llhw    = zeros(1, DIM); 
    
    for k = 1:CLS
        res   = repmat(multi(:,k), 1, DIM) .* x;
        llhw  = llhw - sum(res(y==k,:), 1);
    end
    
    % save the derivative of w_c
    tmpw(:,c) = llhw';
end

% compute the derivative of the objective
tmpw = tmpw + (vertcat(zeros(1,CLS), repmat(reg_w,DIM-1,1)) .* w);
dv   = cat(2, tmpw(:));

 