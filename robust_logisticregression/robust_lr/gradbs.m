% function definition for minFunc
% fv  = function value
% dfv = gradient of the function w.r.t W

function [fv dv] = gradbs(sqb, w, K, y, g, mask, reg_w, reg_b)

% compute regularisation term 
CLS    = size(g, 2);
DIM    = size(K{1}, 2);

w      = reshape(w, DIM, CLS);
reg_w  = reshape(reg_w, 1, CLS);

% compute regularisation term 
reg    = 0.5 * (sum(reg_w .* sum(w(2:end,:))).^2) +  sum(reg_b .* sqb.^2); 

x     = compositeKernel(sqb.^2, K);

% calculate function value
sk     = posteriorYHat(x, w, g);    
fv     = -sum(log(sk(mask)), 1) + reg;

% compute the derivative of the objective


len   =  length(sqb);
dv    =  zeros(len,1);
reg   =  reg_b .* 2 .* sqb;
logit =  softmax(x * w);
for i=1:len
    tmpd = 0;
    for c = 1:CLS        
        % calculate logit w.r.t to w_c and repmat to match sk
        logit_c = repmat(logit(:,c), 1, CLS);
        
        % calculate the big chunk in parentheses
        gc     = repmat(g(c,:), CLS, 1);
        gc     = gc - g;
        chunk  = logit * gc;
        
        % calculate multiplier to X
        tmpw   = sum((logit_c .* chunk) ./ sk, 2) ;
        wt     = w(:,c) * 2 * sqb(i);
        tmpd   = tmpd - sum(tmpw .* (K{i} * wt), 1) + reg(i); 
    end
    dv(i) = tmpd;
end
    







