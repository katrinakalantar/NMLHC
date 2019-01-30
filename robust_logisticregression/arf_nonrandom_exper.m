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

%s = RandStream('mcg16807','Seed',121121); RandStream.setGlobalStream(s);

% global options
ITER            = 10; 
options.maxIter = 100;

EXP_RANGE       = [10 20 30 40 50 100 200 300 400 500 1000];
CLS             = 2;
DIM             = 10;


for nt = {'random'} %'gam12','gam22','gam32','gam51',
    
    noiseType = cell2mat(nt);
    
    
    for i = 1:length(EXP_RANGE)
        
        for j = 1:ITER
            
            [x y ff xx tt dd] = genData(2,DIM,EXP_RANGE(i),1000,2.5,'dis');
                        
            
            Xt = x;
            yt = y;
            Xs = xx;
            ys = tt;
            
            
            [Xt, Xs] = standardise(Xt,Xs);
            
            %% add random/non-random label noise
            target = randi([1 CLS+1])-1;
            
            winit = randn(DIM+1,1);
            
            % train a gold standard decision boundary for noise generation
            % using LR
            options.estG    = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            [ww, ~, ~]      = rlr(winit,eye(CLS),addbias(Xt),yt,options);
            
            
%             %% using SVM
%             ytsvm = castLabel(yt,-1);
%             bestcv = 0;
%             for log2c = -10:10
%                 %for log2g = -10:2:10
%                 cmd = ['-v 5 -c ', num2str(2^log2c)];%, ' -g ', num2str(2^log2g)];
%                 cv = svmtrain(ytsvm, Xt, cmd);
%                 if (cv >= bestcv),
%                     bestcv = cv; bestc = 2^log2c;% bestg = 2^log2g;
%                 end
%                 %end
%             end
%             model = svmtrain(ytsvm, Xt, ['-c ',num2str(bestc)]);%, ' -g ', num2str(bestg)]);
%             w = model.SVs' * model.sv_coef;
%             b = -model.rho;
%             
%             if model.Label(1) == -1
%                 w = -w;
%                 b = -b;
%             end
%             ww = [b; w];
            
            switch noiseType
                case 'random'
                    rate = 0.3;
                    switch target
                        case 1
                            [yz fdz] = injectLabelNoise(yt, [rate 0]);
                        case 2
                            [yz fdz] = injectLabelNoise(yt, [0 rate]);
                        case 0
                            [yz fdz] = injectLabelNoise(yt, [rate rate]);
                    end
                case 'pure'
                    yz = yt;
                    fdz = zeros(size(yz));
                case 'gam12'
                    [yz fdz] = injectNonRandomLabelNoise(yt, ones(size(yt))*-1, addbias(Xt),ww,target,@gammapdf, 1, 2);
                case 'gam22'
                    [yz fdz] = injectNonRandomLabelNoise(yt, ones(size(yt))*-1, addbias(Xt),ww,target,@gammapdf, 2, 2);
                case 'gam32'
                    [yz fdz] = injectNonRandomLabelNoise(yt, ones(size(yt))*-1, addbias(Xt),ww,target,@gammapdf, 3, 2);
                case 'gam51'
                    [yz fdz] = injectNonRandomLabelNoise(yt, ones(size(yt))*-1, addbias(Xt),ww,target,@gammapdf, 5, 1);
                    
            end
            
            % rLR
            options.estG = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            [w g l_lr] = rlr(winit,eye(2),addbias(Xt),yz,options);
            es_lr(j,i) = sum(sign(addbias(Xs)*w) ~= castLabel(ys,-1))/length(ys);
            
            % rLR
            options.estG = true;
            options.regFunc = 'lasso';
            options.verbose = false;
            rr = .2;
            [wr gr l_rlr] = rlr(winit,[1-rr rr;rr 1-rr],addbias(Xt),yz,options);
            es_rlr(j,i) = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
            
            % nLR
            options.estG = true;
            options.regFunc = 'lasso';
            options.verbose = false;
            [wg nd l_gam] = gammalr(winit,addbias(Xt), yz, options);
            es_nlr(j,i) =  sum(sign(addbias(Xs)*wg) ~= castLabel(ys,-1))/length(ys);
            
            
            
            fprintf('[%d %d] lr=%f rlr=%f nlr=%f\n',i,j,es_lr(j,i),es_rlr(j,i),es_nlr(j,i));
            
        end
                
    end
    
    disp('============================================================');
    % some data information
    fprintf(1,'Stat: varied dataset size\n');
    for i=1:size(es_rlr,2)
        fprintf('[%d %d] rlr=%f gam=%f hlr=%f\n',i,j,...
            mean(es_lr(:,i)), mean(es_rlr(:,i)), mean(es_nlr(:,i)));
    end
    disp('------------------------------------------------------------');
    for i=1:size(es_rlr,2)
        fprintf('[%d %d] rlr=%f gam=%f hlr=%f\n',i,j,...
            median(es_lr(:,i)), median(es_rlr(:,i)), median(es_nlr(:,i)));
    end
    disp('============================================================');
    
    %graphFactory(ITER, EXP_RANGE, '','Dataset Size','Generalisation Error (%)',...
    %    es_lr.*100,'LR',...
    %    es_rlr.*100,'rLR',...
    %    es_nlr.*100,'gLR');
    
    errorbar(EXP_RANGE,mean(es_lr.*100),std(es_lr.*100)); 
    errorbar(EXP_RANGE,mean(es_rlr.*100),std(es_rlr.*100)); 
    errorbar(EXP_RANGE,mean(es_nlr.*100),std(es_nlr.*100));
    xlabel('Dataset Size');
    ylabel('Generalisation Error (%)');
    
    %axis tight;
    ylim([0,50]);
    %set(gca,'xtickmode','manual','yminortick','off');
    set(gca,'ytick',[0:10:50]);
    % extra for log scale
    set(gca,'xscale','log');
    
    %if strcmp(nt,'random');
    %    legend('Location','NorthWest');
    %else
    legend('Location','NorthEast');
        
    %end
    filename = strcat('/Users/kkalantar/Downloads/',nt,'_vary_data');
    %fprintf(filename);
    figname  = strcat(filename,'.eps');
    print(gcf, '-dpsc2', figname{1});
    
end


