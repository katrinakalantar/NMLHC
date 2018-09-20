% Name  :  Robust Multinomial Logistic Regression
% Author:  Jakramate Bootkrajang
% Last update: 16 April 2013
% Input :  
%          w       = given initial model's parameters  (DIM,CLS)
%          g       = initial label flipping probability matrix
%          x       = Design matrix where a row represents a sample
%          y       = Target value
%       options.maxIter = maximum optimising iteration
%       options.estGam  = estimating the gamma using multiplicative update
%       options.regFunc = type of regularisation: 'noreg', 'L1' or 'L2'
%       options.sn      = Small number for approximating L1 objective.
%       options.verbose = displaying negative log-likelihood
% Output:  
%          w     = fitted model
%          g     = estimated gamma matrix
%          llh   = negative log-likelihood
% Note  :
%       1. The model uses {1,2,...K} class representation internally. 
%  |-----|-----------|
%  |     |  1   y^  n|
%  ------|-----------|
%  |y  1 | g00   g0n |   where n = k-1 and k = number of class
%  |   n | gn0   gnn |
%  |-----|-----------|
%=====================================================================================   

function [w g llh] = rmlr(w, g, x, y, options)

% check if bias term is added i.e. first 
% column must be unique and equal to 1
if (sum(x(:,1)) ~= size(x,1))
    disp('Bias terms might have not been added');
end

if ~isfield(options,'sn')
    options.sn = 1e-8;
end

if ~isfield(options,'maxIter')
    options.maxIter = 50;
end

CLS    = length(unique(y)); 
DIM    = size(x,2);

% create label mask 
mask = repmat(y, 1, CLS);
for k = 1:CLS
    mask(:,k) = logical(y==k);
end
mask = logical(mask);

% ========================= BEGIN ESTIMATING PARAMETERS =========================== 
for l=1:options.maxIter


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
    
    switch options.regFunc
        case 'lasso'
            lambda      = 1 ./(sqrt(w.^2 + options.sn));    
            lambda(1,:) = 0;             
        case 'l2'                       
            lambda      = (length(w)./2 +1)./((sum(sum(w(2:end,:))).^2)./2 + 2); 
            lambda(1,:) = 0; 
        case 'noreg'
            lambda      = zeros(size(w));
        otherwise
            disp('Invalid regularisation function. Now using L1');
            options.regFunc = 'L1';
    end
    if l==1
        lambda = lambda * 0;
    end   
    
    % minimise objective function 
    w0               = cat(2, w(:));
    lambda0          = cat(2, lambda(:));

    opts.maxIter     = 10;
    opts.Display     = false;    
    [w0, fv, v1, v2] = minFunc(@gradws, w0, opts, g, x, y, mask, lambda0, options);               
    w                = reshape(w0, DIM, CLS);             
    
    % save the progress
    llh(l)           = fv(end);
    
    if options.verbose
        disp(llh(l));
    end    
    
end % end optimisation loop    
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


