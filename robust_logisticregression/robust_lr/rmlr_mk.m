% Name  :  Robust Multi-Kernel Multinomial Logistic Regression
% Author:  Jakramate Bootkrajang
% Last update: 16 April 2013
% Input :
%          w    = initial weight vector
%          g    = initial gamma matrix
%          b    = initial kernel coefficient
%          K    = Cell of kernels
%          y    = Target value
%          options.maxIter = max iterations
%          options.estGam  = estimate gamma ? {true,false}
%          options.estBeta = estimate beta  ? {true,false}
%          options.reg     = add regularisation ? {default=true}
%          options.verbose = {true,false}
% Output:
%          w    = fitted model
%          g    = estimated gamma matrix
%          b    = estimated beta
%          llh  = model's log likelihood
% Note  :
%       1. The model uses {1,2,...K} class representation internally.
%=====================================================================================   

function [w g b llh] = rmlr_mk(w, g, b, K, y, options)

if ~isfield(options,'maxIter')
    options.maxIter = 50;
end

CLS    = size(unique(y),1); 
DIM    = size(K{1},2);

% create label mask 
mask = repmat(y, 1, CLS);
for k = 1:CLS
    mask(:,k) = logical(y==k);
end
mask = logical(mask);


% ========================= BEGIN ESTIMATING PARAMETERS =========================== 
for i=1:options.maxIter 

    % calculate regularisation 
    if options.reg
        reg_ws  = (length(w)./2 +1)./((sum(w(2:end,:)).^2)./2 + 2);         
        reg_b   = (1 + 1)./(b + 1e-100); 
    else
        reg_ws  = 0;        
        reg_b   = 0;
    end
    
    x = compositeKernel(b, K);
    
    if options.reg
        opts.MaxIter = 10;
    else
        opts.MaxIter = 100;
    end    
    % minimise objective function 
    w0               = cat(2, w(:));
    reg_ws0          = cat(2, reg_ws(:));     
    [w0, fv, v1, v2] = minFunc(@gradws_mk, w0, opts, g, x, y, mask, b, reg_ws0, reg_b);               
    w                = reshape(w0, DIM, CLS);   
    
    if i>1;  llh(i) = fv(end); end;
    
    if options.reg
        reg_ws  = (length(w)./2 +1)./((sum(w(2:end,:)).^2)./2 + 2); 
    else
        reg_ws  = 0;
    end    
    
    if options.estB
        opts.MaxIter = 3;
        [sqb, fb, vb, eb] = minFunc(@gradbs,sqrt(b),opts,w,K,y,g,mask,reg_ws,reg_b);  
        b   = sqb.^2;  
    end
    
    x = compositeKernel(b, K);
                
    if options.estG
        for gl=1:3
            % update the gamma
            sk = posteriorYHat(x, w, g);
            lj = softmax(x * w);
            num = zeros(CLS,CLS);
            for k = 1:CLS
                skc = repmat(sk(:,k), 1, CLS);
                tmp = lj ./ skc;
                % select y_hat_i = k with mask
                tmp2 = tmp(mask(:,k),:);
                num(:,k) = sum(tmp2,1);
            end
            
            % calculate denominator
            dnt   = sum((num .* g), 2);
            denom = repmat(dnt, 1, CLS);
            
            % now last step for gamma
            g     = g .* (num ./ denom);
        end
    end                  

end % end optimisation loop    

llh = llh(2:end);
% ============================== END ESTIMATING PARAMETER ============================

% class identity assignment, solved with Hungarian algorithm
task    = repmat(max(max(g)),CLS,CLS) - g;
[asm n] = munkres(task);
s       = 1;

for loop = 1:CLS
    [nu i] = max(asm(s,:));
    if (asm(s,i) > asm(i,i))
        asm = rowswap(asm, s, i);
        g   = rowswap(g, s, i);        
        w   = rowswap(w',s, i)';
    else
        s   = s + 1;
    end
end


