% =========================================================================
% A wrapper experiment script for the robust Logistic Regression algorithm
% =========================================================================

addpath('./utils/');
addpath('./kk_utils/');
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

%addpath('/Users/kkalantar/Documents/Downloads/PCAMV');
 
% start with nice & clean workspace
clear all;
close all;

% global options
ITER            = 10; %10;%10; 
options.maxIter = 100;

% parameters iterated
EXP_RANGE       = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5]; %y-axis of heatmap
EXP_RANGE_J     = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap
%EXP_RANGE       = [10 20 30 40 50 100 200 300 400 500 1000];%this range was used for sample sizes
%EXP_RANGE       = [10 50 100 500 1000 2500 5000];% 10000 20000];%this range was used for sample sizes
%EXP_RANGE_J     = [.3];%, .05, .1, .15, .2, .25, .3];%, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap

CLS             = 2;    % number of classes
DIM             = 100; % dimensionality of dataset generated
DS_SIZE         = 100;  % dataset size

% preallocating error storage
es_lr      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_lr_nonoise      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_gammalr = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);

% preallocating AUC storage
es_lr_auc       = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr_auc      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_gammalr_auc  = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_lr_nonoise_auc = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);

index = 0;
use_PCs = false;
feature_select = false; %true;
%n_features = 2000;


for i = 1:length(EXP_RANGE)
    
    flip_i = EXP_RANGE(i);
    %flip_i = 0;
    %DIM = EXP_RANGE(i); %10000;
    n_features = DIM;
    %DS_SIZE = 100;% EXP_RANGE(i);
    
    for j = 1:length(EXP_RANGE_J)
        
        flip_j = EXP_RANGE_J(j);
        index = index + 1;
        
        for k = 1:ITER
            
            [x, y, ff, xx, tt, dd] = genData(2,DIM,DS_SIZE,1000,2.5,'gen'); %second "DIM" used to be "1000"
            %[x y fd xx tt dd] = genData(2,5,500,10,0.5,'gen');
            %[x, y, ff, xx, tt, dd] = subsetData(dataset_x, dataset_y, .8);

            Xt = double(x);
            yt = y;
            Xs = double(xx);
            ys = tt;
            
            [Xt, Xs] = standardise(Xt,Xs);  % get a "step direction is illegal" error if you don't standardize
            
            size(Xt)
            
            [yz, fdz] = injectLabelNoise(yt, [flip_i flip_j]);
            fprintf('[i = %f j = %f] sum=%f flip_ids=%f\n',flip_i,flip_j,sum(fdz ~= -1),(flip_i*length(yt) + flip_j*length(yt))/2);

            % add random/non-random label noise
            target = randi([1 CLS+1])-1;

            % create a new winit for training the noised, feature-selected !
            winit = randn(size(Xt,2) + 1, 1);
            
            common_reg = 'lasso';
            common_sn = 1e-8;
            
            
            % rLR (using true labels) to give baseline of number of
            % features
%             options.estG = false;
%             options.regFunc = common_reg;
%             options.verbose = false;
%             options.sn = common_sn;
%             [w, g, l_lr] = rlr(winit,eye(2),addbias(Xt),yt,options);
%             [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*w,2);
%             es_lr_auc(i,j,k) = AUC;
%             sum(abs(w) < 1e-6)/length(w)
            
            % rLR
            options.estG = false;
            options.regFunc = common_reg;
            options.verbose = false;
            options.sn = common_sn;
            [w, g, l_lr] = rlr(winit,eye(2),addbias(Xt),yz,options);
            es_lr(i,j,k) = sum(sign(addbias(Xs)*w) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*w,2);
            es_lr_auc(i,j,k) = AUC;
            
            %sum(abs(w) < 1e-6)/length(w)
            
            % rLR
            options.estG = true;
            options.regFunc = common_reg;
            options.verbose = false;
            options.sn = common_sn;
            rr = .2;
            [wr, gr, l_rlr] = rlr(winit,[1-rr rr;rr 1-rr],addbias(Xt),yz,options);
            es_rlr(i,j,k) = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr,2);
            es_rlr_auc(i,j,k) = AUC;
            
            %sum(abs(wr) < 1e-6)/length(wr)
           
            % nLR
            options.estG = true;
            options.regFunc = common_reg;
            options.verbose = false;
            options.sn = common_sn;
            [wg, nd, l_gam] = gammalr(winit,addbias(Xt), yz, options);
            es_gammalr(i,j,k) =  sum(sign(addbias(Xs)*wg) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wg,2);
            es_gammalr_auc(i,j,k) = AUC;
            
            %sum(abs(wg) < 1e-6)/length(wg)
            
            %fprintf('AUC lr=%f rlr1=%f rlr2=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_rlr_auc2(i,j,k),es_gammalr_auc(i,j,k));
            fprintf('AUC lr=%f rlr1=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_gammalr_auc(i,j,k));
            fprintf('ERR lr=%f rlr1=%f gammalr=%f\n',es_lr(i,j,k),es_rlr(i,j,k),es_gammalr(i,j,k));
            
            
        end
            
    end
   
end

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure6','-bestfit','-dpdf')

plot_heatmap(es_lr, 'MEAN ES\_LR', 'figure1', 0, .55, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_rlr, 'MEAN ES\_RLR', 'figure2', 0, .55, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr, 'MEAN ES\_GAMMALR', 'figure3', 0, .55, EXP_RANGE_J, EXP_RANGE)

plot_heatmap(es_lr_auc, 'MEAN ES\_LR AUC', 'figure1_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_rlr_auc, 'MEAN ES\_RLR AUC', 'figure2_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr_auc, 'MEAN ES\_GAMMALR AUC', 'figure3_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)

plot_lineplot(es_lr, es_rlr, es_gammalr ,EXP_RANGE_J, EXP_RANGE, 'Figure4.pdf')
%plot_lineplot(100 - mean(es_lr,3), 100 - mean(es_rlr,3), 100 - mean(es_gammalr,3) ,EXP_RANGE_J, EXP_RANGE, 'Figure5')
plot_lineplot((1 - es_lr), (1 - es_rlr), (1 - es_gammalr) ,EXP_RANGE_J, EXP_RANGE, 'Figure5')
plot_lineplot(es_lr_auc, es_rlr_auc, es_gammalr_auc ,EXP_RANGE_J, EXP_RANGE, 'Figure5_auc')


for i = 1:size(es_lr,1)
    p = ranksum(reshape(es_lr(i,:,:),[],1),reshape(es_rlr(i,:,:),[],1));
    final_val = "";
    if p < .05
        final_val = "*";
    end
    fprintf("p = %f at value = %f   %s\n", p, EXP_RANGE(i), final_val)
end

for i = 1:size(es_lr,1)
    p = ranksum(reshape(es_lr(i,:,:),[],1),reshape(es_gammalr(i,:,:),[],1));
    final_val = "";
    if p < .05
        final_val = "*";
    end
    fprintf("p = %f at value = %f   %s\n", p, EXP_RANGE(i), final_val)
end
