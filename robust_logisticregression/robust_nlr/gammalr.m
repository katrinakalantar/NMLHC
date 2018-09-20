%==========================================================================
% Name  :  Robust Logistic Regression + Gamma Noise function
% Author:  Jakramate Bootkrajang
% Last update: 20 February 2016
% Input :
%       w = given initial model's parameters
%       x = Design matrix where a row represents a sample
%       y = Target value
%       options.maxIter = maximum optimising iteration
%       options.regFunc = type of regularisation: 'noreg', 'L1' or 'L2'
%       options.sn      = Small number for approximating L1 objective.
%       options.verbose = displaying negative log-likelihood
%       options.estND   = estimating noise distribution
% Output:
%       w   = fitted model
%       nd  = estimated noise distribution
%       llh = negative log-likelihood
% Note  :
%       1. The function uses {0,1} class representation.
%       w = (eye(size(x,2))-x'*inv(x*x'+eye(size(x,1)))*x)*x'*y/2;  w(1) = 1;
%  |-----|-----------|
%  |     |  0   y^  1|
%  ------|-----------|
%  |y  0 | g00   g01 |
%  |   1 | g10   g00 |
%  |-----|-----------|
%==========================================================================
%==========================================================================

function [w nd llh] = gammalr(w, x, y, options)


% The function uses {0,1} class representation.
y = castLabel(y,0);

%% Some input preprocessing
% check if the problem is binary problem
if (size(unique(y),1) ~= 2) && options.verbose
    error('Trying to solve multiclass problem');
end

% checking if bias term is added i.e. first
% column must be unique and equals to 1
if (sum(x(:,1)) ~= size(x,1)) && options.verbose
    disp('Bias terms might have not been added');
end

if ~isfield(options,'maxIter')
    options.maxIter = 50;
end

if ~isfield(options,'sn')
    options.sn = 1e-8;
end

if ~isfield(options,'estND')
    options.estND = true;
end


% parameters to the gamma function
t0 = 500;
t1 = 500;
k0 = 1;
k1 = 1;

% for storing log-likelihood values
llh = zeros(1,options.maxIter);
l   = 1;

optsw.maxIter = 5;
optsw.Display = false;
optso.maxIter = 3;
optso.Display = false;


%% ========================= BEGIN ESTIMATING PARAMETERS ===========================
while l<=options.maxIter
        
%     options.bReg0 = 2/(log(t0-1)+10^8);%3*(log(t0-1)-mu)/((log(t0-1)-mu)^2 + 2*beta);
%     options.bReg1 = 2/(log(t1-1)+10^8);%3*(log(t1-1)-mu)/((log(t1-1)-mu)^2 + 2*beta);

    options.bReg0 = -1/(1-t0);
    options.bReg1 = -1/(1-t1);
    
    % no regularisation and use uniform noise on the first iter
    if l==1
        lambda = zeros(size(w));        
        g01    = ones(size(y)) * 0.2;
        g10    = ones(size(y)) * 0.2;
    else        
        g01    = gammapdf(z,1,t0);
        g10    = gammapdf(z,1,t1);
    end
    
    
    % optimising the weight vector         
    [w, ~, ~, ~] = minFunc(@gradw_gammalr, w, optsw, g01, g10, x, y, t0, t1, lambda, options);            
    %[w, ~, ~, ~] = minFunc(@gradw_nr_new, w, optsw, x, y, t0, t1, lambda, options);            
    %dw = checkgrad('gradw_nr_new',w, 1e-7, x, y, t0,t1,lambda, options)    
         
    
    
    % update z based on new weight vector
    t = x * w;     
    if l < 3 
        z = t/norm(w);
        %z = corr(x',w);  
        
    end
    
    
    %% recomputing the regularisation term based on new w
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
    
    [regV, ~]    = feval(options.regFunc, w, options.sn);
    options.regV = sum(lambda .* regV);
         
    
    % estimating noise parameters
    if options.estG
        
        %dt0 = checkgrad('gradt0',t0, 1e-7, z, t, y, t1, k0, k1, options)
        %dt1 = checkgrad('gradt1',t1, 1e-7, z, t, y, t0, k0, k1, options)
        
        [t0,  ~, ~, ~] = minFunc(@gradt0, t0, optso, z, t, y, t1, k0, k1, options);
        [t1, fv, ~, ~] = minFunc(@gradt1, t1, optso, z, t, y, t0, k0, k1, options);                                            
       
    end

    % calculating likelihood function           
    llh(l) = fv;         
    l  =  l + 1; 
 
    % if you need extra information
    if options.verbose                
        subplot(4,1,4);
        plot(llh(2:l-1),'bx','LineWidth',2); drawnow;
        axis tight;
    end

end

%% ======================== END ESTIMATING PARAMETER =======================

nd.t0 = t0;
nd.t1 = t1;
nd.k0 = k0;
nd.k1 = k1;

