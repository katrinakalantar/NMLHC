% =========================================================================
% A wrapper experiment script for the robust Logistic Regression algorithm
% =========================================================================

addpath('./utils/');
addpath('./robust_nlr/');
addpath('./robust_lr/');
addpath('./regulariser/');
addpath('./third_party_libs/');
addpath('./third_party_libs/minFunc');
addpath('./third_party_libs/netlab');
addpath('./third_party_libs/libsvm');
addpath('./third_party_libs/glmnet_matlab');
addpath('./third_party_libs/minFunc/mex');
addpath('./third_party_libs/minFunc/compiled');

% start with nice & clean workspace
clear all
close all

%s = RandStream('mcg16807','Seed',121121); RandStream.setDefaultStream(s);

% global options
ITER            = 10;%0;
options.maxIter = 100;

% type of noises
nt              = {'pure','gam12','gam22','gam32','gam51','random'};
EXP_RANGE       = 1:length(nt);


% a grid for visualisation
gridRange = -20:0.05:20;
posGrid   = 0:0.05:5;
negGrid   = -5:0.05:0;
[x1,x2]   = meshgrid(gridRange, gridRange);
xgrid     = [reshape(x1,size(x1,1)^2,1) reshape(x2,size(x2,1)^2,1)];

%balance data
%for dat = {}
%imbalance data
for dat = {'dataset/colon','dataset/leukaemia','dataset/websearch','dataset/breast'}
    
    % load real data
    data = cell2mat(dat); load(data);
    
    % generate artificial data
    %[x y fd xx tt dd] = genData(2,5,500,10,0.5,'gen');
    
    
    % fixing the labels if ground truth is known
    yori = y;                                 % yori are the original labels
    y    = castLabel(correctLabel(y,fd),0);   % y are the corrected labels
    CLS  = length(unique(y));
    
    % probing data dimension and some info
    [DPnt DIM] = size(x);
    
    
    % set training set size
    if strcmp(data,'breast');
        TrPnt  = floor(0.9*DPnt);
    else
        TrPnt  = floor(0.8*DPnt);
    end
    
    % preallocating error storage
    es_lr      = nan(ITER,max(EXP_RANGE));
    es_rlr     = nan(ITER,max(EXP_RANGE));
    es_gammalr = nan(ITER,max(EXP_RANGE));
    
    clsf = who('es_*');
    
    
    
    for i = EXP_RANGE
        % looping through each noise type
        noiseType = cell2mat(nt(i));
        
        for j = 1:ITER
            % random permutation of the dataset
            perm = randperm(DPnt);
            if strcmp(noiseType,'pure')
                Xt = x(perm(1:TrPnt),:);
                yt = yori(perm(1:TrPnt));
                Xs = x(perm(TrPnt+1:end),:);
                ys = y(perm(TrPnt+1:end));
            else
                Xt = x(perm(1:TrPnt),:);
                yt = y(perm(1:TrPnt));
                Xs = x(perm(TrPnt+1:end),:);
                ys = y(perm(TrPnt+1:end));
            end
            
            % normalising the data
            [Xt, Xs] = standardise(Xt,Xs);
            
            % adding random/non-random label noise
            target = randi([1 CLS]);
            
            
            %             % training a gold standard decision boundary for noise generation
            %             % using LR
            winit = zeros(DIM+1,1);
            options.estG    = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            [wstar, ~, ~]   = rlr(winit,eye(CLS),addbias(Xt),yt,options);
            
            
%                         %% using SVM
%                         ytsvm = castLabel(yt,-1);
%                         bestcv = 0;
%                         for log2c = -10:10
%                             for log2g = -10:2:10
%                             cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
%                             cv = svmtrain(ytsvm, Xt, cmd);
%                             if (cv >= bestcv),
%                                 bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
%                             end
%                             end
%                         end
%                         model = svmtrain(ytsvm, Xt, ['-c ',num2str(bestc), ' -g ', num2str(bestg)]);
%                         w = model.SVs' * model.sv_coef;
%                         b = -model.rho;
%             
%                         if model.Label(1) == -1
%                             w = -w;
%                             b = -b;
%                         end
%                         wstar = [b; w];
            
            
            if strcmp(noiseType,'random')
                rate = 0.3;
                switch target
                    case 1
                        [yz fdz] = injectLabelNoise(yt, [rate 0]);
                    case 2
                        [yz fdz] = injectLabelNoise(yt, [0 rate]);
                    case 0
                        [yz fdz] = injectLabelNoise(yt, [rate rate]);
                end
                shape = 1; scale=1;
            elseif strcmp(noiseType, 'pure')
                yz = yt;
                fdz = zeros(size(yz));
            else
                switch noiseType
                    case 'gam12'
                        shape = 1; scale=2;
                        func = 'gammapdf';
                    case 'gam22'
                        shape = 2; scale=2;
                        func = 'gammapdf';
                    case 'gam32'
                        shape = 3; scale=2;
                        func = 'gammapdf';
                    case 'gam51'
                        shape = 5; scale=1;
                        func = 'gammapdf';
                end
                [yz fdz] = injectNonRandomLabelNoise(yt, ones(size(yt))*-1, addbias(Xt), wstar, target,str2func(func), shape, scale);
            end
            
            
            % LR-benchmark
            es_lr(j,i) = sum(sign(addbias(Xs)*wstar) ~= castLabel(ys,-1))/length(ys);
            
            commonReg = 'lasso';
            
            
            % rLR
            options.estG      = true;
            options.regFunc   = commonReg;
            options.verbose   = false;
            gamMat            = [0.8 0.2;0.2 0.8];
            [wr gamMat l_rlr] = rlr(winit,gamMat,addbias(Xt),yz,options);
            es_rlr(j,i)       = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
            
            % gammaLR
            options.estG    = true;
            options.regFunc = commonReg;
            options.verbose = false;
            [wga nd l_gam]  = gammalr(winit, addbias(Xt), yz, options);
            es_gammalr(j,i) = sum(sign(addbias(Xs)*wga) ~= castLabel(ys,-1))/length(ys);
            
            
            fprintf('[%d %2d] ',i,j);
            for c=1:length(clsf)
                curStr = cell2mat(clsf(c));
                curE   = eval(curStr);
                fprintf('%s=%2.3f ',curStr(4:end), curE(j,i));
            end
            if target==0
                fprintf('**\n');
            else
                fprintf('\n');
            end
            
        end
        
    end
    
    disp('==================================================================');
    % some data information
    fprintf(1,'Stat: data = %s, cardinality = %d balanceness = %2.2f\n',data, size(x,2), sum(y==1)/length(y));
    disp('----------------------------++ MEAN ++----------------------------');
    for i=1:max(EXP_RANGE)
        fprintf('[%d %2d] ',i,j);
        for c=1:length(clsf)
            curStr = cell2mat(clsf(c));
            curE   = eval(curStr);
            fprintf('%s=%2.5f ',curStr(4:end), mean(curE(:,i)));
        end
        fprintf('\n');
    end
    disp('---------------------------++ MEDIAN ++---------------------------');
    for i=1:max(EXP_RANGE)
        fprintf('[%d %2d] ',i,j);
        for c=1:length(clsf)
            curStr = cell2mat(clsf(c));
            curE   = eval(curStr);
            fprintf('%s=%2.5f ',curStr(4:end), median(curE(:,i)));
        end
        fprintf('\n');
    end
    disp('==================================================================');
    
    
end

