% Name  :  Robust Multi-Kernel Logistic Regression
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
%       1. The model uses {0,1} class representation internally.
%==========================================================================

function [w g b llh] = rlr_mk(w, g, b, K, y, options)

% check if the problem is binary problem
if (size(unique(y),1) ~= 2)
    error('Trying to solve multiclass problem');
end

if ~isfield(options,'maxIter')
    options.maxIter = 50;
end

% suppress minFunc's output
opts.Display = 'off';

% =================== BEGIN ESTIMATING PARAMETERS =========================
for i=1:options.maxIter
    
    % calculate regularisation
    if options.reg
        reg_w  = (length(w)/2 +1)/(sum(w(2:end).^2)/2 + 2);
        reg_b  = (1 + 1)./(b + 1e-100);
    else
        reg_w  = 0;
        reg_b  = 0;
    end
    
    x = compositeKernel(b, K);
    
    if options.reg
        opts.MaxIter = 10;
    else
        opts.MaxIter = 100;
    end
    [w, fw, vw, ew] = minFunc(@gradw_mk,w,opts,x,y,g,b,reg_w,reg_b);
    
    %checkgrad('gradw_mk', w, 1e-5, x,y,g,b,reg_w,reg_b)
    
    if i>1;  llh(i) = fw(end); end;
    
    if options.reg
        reg_w  = (length(w)/2 + 1)/(sum(w(2:end).^2)/2 + 2);
    else
        reg_w  = options.ureg;  % user supplied regularisation parameter
    end
    
    if options.estB
        opts.MaxIter = 3;        
        [sqb, fb, vb, eb] = minFunc(@gradb,sqrt(b),opts,w,K,y,g,reg_w,reg_b);
        b = sqb.^2;
    end
    
    x = compositeKernel(b, K);
    t = x * w;
    
    if options.estG && i < 4
        for gl=1:1
            %for gl=1:1
            s0  = (g(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,1) ./ (1+exp(-t)));
            s1  = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));
            
            % avoid numerical problem
            s0(s0==0) = eps;   
            s1(s1==0) = eps;
            
            % update the gammas
            g00 = g(1,1) * sum(((1 - y) ./ (s0)) .* (1 ./ ((1 ./ exp(-t))+1)),1);
            g01 = g(1,2) * sum((y ./ s1) .* (1 ./ ((1 ./ exp(-t))+1)),1);
            
            % (g00 + g01) is the Lagrangian
            g(1,1) = g00 / (g00 + g01);
            g(1,2) = 1 - g(1,1);
            
            g10 = g(2,1) * sum(((1-y) ./ s0) ./(1+exp(-t)),1);
            g11 = g(2,2) * sum(((y) ./ (s1)) ./(1+exp(-t)),1);
            
            % (g10 + g11) is the Lagrangian
            g(2,1) = g10 / (g10 + g11); 
            g(2,2) = 1 - g(2,1);
        end
    end
end

llh = llh(end);

% ===================== END ESTIMATING PARAMETER ==========================

% deal with class identity problem
s = 1;
for i = 1:2
    [vm, m] = max(g(s,:));
    if (g(s,m) > g(m,m))
        g = circshift(g,1);
        w = -w;
    else
        s = s + 1;
    end
end
