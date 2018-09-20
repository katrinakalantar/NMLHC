%==========================================================================
% Name  :  Robust Logistic Regression
% Author:  Jakramate Bootkrajang
% Last update: 16 April 2013
% Input :
%       w = given initial model's parameters
%       g = initial label flipping probability matrix
%       x = Design matrix where a row represents a sample
%       y = Target value
%       options.maxIter = maximum optimising iteration
%       options.estGam  = estimating the gamma using multiplicative update
%       options.regFunc = type of regularisation: 'noreg', 'L1' or 'L2'
%       options.boost   = adapting to boosting framework
%       options.sn      = Small number for approximating L1 objective.
%       options.verbose = displaying negative log-likelihood
% Output:
%       w   = fitted model
%       g   = estimated gamma matrix
%       llh = negative log-likelihood
% Note  :
%       1. The function uses {0,1} class representation.
%       w = (eye(size(x,2))-x'*inv(x*x'+eye(size(x,1)))*x)*x'*y/2;  w(1) = 1;
%  |-----|-----------|
%  |     |  0   y^  1|
%  ------|-----------|
%  |y  0 | g00   g01 |
%  |   1 | g10   g11 |
%  |-----|-----------|
%==========================================================================
%==========================================================================

function [w g llh] = rlr(w, g, x, y, options)


% The function uses {0,1} class representation.
y = castLabel(y,0);

%% Some input preprocessing
% check if the problem is binary problem
if (size(unique(y),1) ~= 2) && options.verbose
    error('Trying to solve multiclass problem');
end

% checking if bias term is added i.e. first
% column must be unique and equal to 1
if (sum(x(:,1)) ~= size(x,1)) && options.verbose
    disp('Bias terms might have not been added');
end

if ~isfield(options,'sn')
    options.sn = 1e-8;
end

if ~isfield(options,'maxIter')
    options.maxIter = 50;
end

if ~isfield(options,'boost')
    options.boost = false;
end

if ~isfield(options,'estG')
    options.estG = true;
end

if ~isfield(options,'optimType')
    options.optimType = 1;
end

if options.boost
    D = options.dist;
else
    D = ones(size(x,1),1);
    options.dist = D;
end

% for storing log-likelihood values
llh = zeros(1,options.maxIter);


%% ========================= BEGIN ESTIMATING PARAMETERS ===========================
for l=1:options.maxIter
    
    %% computing the regularisation term
    switch options.regFunc
        case 'lasso'
            % LASSO-like regularisation on top of each variable
            % first variable is a bias term, no reg. needed
            lambda    = 1./(sqrt(w.^2 + options.sn));
            lambda(1) = 0;
        case 'l2'
            % Bayesian L2 regularisation
            lambda    = (length(w)/2 +1)/(sum(w(2:end).^2)/2 + 2);
        case 'noreg'
            % No regularisation imposed
            lambda    = 0;
        case 'fixed'
            lambda    = ones(size(w));
            lambda(1) = 0;
        otherwise
            disp('Invalid regularisation function. Now using L1');
            options.regFunc = 'lasso';
    end
    
    % no regularisation in the first run
    if l==1
        lambda = ones(size(w));
        lambda(1) = 0;                    
    end    
    
    t  = x * w;
    
    %% updating the gamma matrix
    if options.estG %&& l < 5
        for gg=1:1
        % using multiplicative update
        s0  = (g(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,1) ./ (1+exp(-t)));
        s1  = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));
        
        % avoiding numerical problem
        s0(s0==0) = eps;
        s1(s1==0) = eps;
        
        g00 = g(1,1) * sum((D.*(1 - y) ./ (s0)) .* (1 ./ ((1 ./ exp(-t))+1)),1);
        g01 = g(1,2) * sum((D.*y ./ s1) .* (1 ./ ((1 ./ exp(-t))+1)),1);
        
        % (g00 + g01) is the Lagrangian
        g(1,1) = g00 / (g00 + g01);
        g(1,2) = 1 - g(1,1);
        
        g10 = g(2,1) * sum((D.*(1-y) ./ s0) ./(1+exp(-t)),1);
        g11 = g(2,2) * sum((D.*(y) ./ (s1)) ./(1+exp(-t)),1);
        
        % (g10 + g11) is the Lagrangian
        g(2,1) = g10 / (g10 + g11);
        g(2,2) = 1 - g(2,1);
        end
    end
       
    
    %% estimating the weight vector
    switch options.optimType
        case 1
            opts.maxIter    = 10;
            opts.Display    = false;
            [w, fv, v1, v2] = minFunc(@gradw, w, opts, g, x, y, lambda, options);
            
        case 2
            % requires rp_eda.m presented in GECCO 2013 paper
            % 'Towards large scale EDA'            
            opts.maxFE    = 200000;
            opts.optmType = 'min';
            opts.ub       = 100;
            opts.lb       = -100;
            opts.topSize  = 75;
            opts.covType  = 'full';
            opts.k        = 3;
            opts.rpmSize  = 1000;
            opts.verbose  = false;
            opts.typeR    = 'G';
            opts.popSize  = 100;
                        
            [w, fv, st, op]  = rp_eda(@rlrobj, opts,...
                randn(opts.popSize,size(w,1)),...
                g, x, y, repmat(lambda,1,opts.popSize), options);
            w = w';
        otherwise
            options.optimType = 1;
    end
    
    if options.verbose        
        if l==1; figure(2); end
        subplot(2,1,1);
        plot(llh(2:l),'bx','LineWidth',2); drawnow;
        
        subplot(2,1,2);        
        tidy(w,'r','xx','yy'); drawnow;
    end
    
    t  = x * w;
    
    s0    = (g(1,1) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,1) ./ (1+exp(-t)));
    s1    = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));
    
    s0(s0==0) = eps;
    s1(s1==0) = eps;
    
    [regV void] = feval(options.regFunc, w, options.sn);
    
    llh(l) = -sum((y .* log(s1)) + ((1 - y).* log(s0)),1) + sum(lambda .* regV);
            
end

if options.verbose
    %close;
end
% ======================== END ESTIMATING PARAMETER =======================

%% dealing with class identity problem
s = 1;
for loop = 1:2
    [un, i] = max(g(s,:));
    if (g(s,i) > g(i,i))
        g = circshift(g,1);
        w = -w;
    else
        s = s + 1;
    end
end

%llh = llh(end);  % return only the last value of llh

