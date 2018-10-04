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
ITER            = 2; %10;%10; 
options.maxIter = 100;

% parameters iterated
%EXP_RANGE       = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5]; %y-axis of heatmap
%EXP_RANGE_J     = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap
EXP_RANGE       = [0, .1, .2,  .3];%, .05, .1, .15, .2, .3, .4, .5]; %y-axis of heatmap.025, 
EXP_RANGE_J     = [0];%, .05, .1, .15, .2, .25, .3];%, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap

CLS             = 2;    % number of classes
DIM             = 1000; % dimensionality of dataset generated
DS_SIZE         = 50;  % dataset size

% preallocating error storage
es_lr      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_lr_nonoise      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr2    = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
% es_rlr3    = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_gammalr = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);

% preallocating AUC storage
es_lr_auc       = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr_auc      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr_auc2     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
% es_rlr_auc3     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_gammalr_auc  = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_lr_nonoise_auc = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);

index = 0;
use_PCs = false;
feature_select = true;
n_features = 2000;

% %%% LOAD BIOLOGY DATA
% gseData = geoseriesread('/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE60244_series_matrix.txt');
% sampleSources = unique(gseData.Header.Samples.source_name_ch1);
% bacteria = gseData.Data(:,contains(string(gseData.Header.Samples.source_name_ch1), 'BACTERIA'));  % bacterial samples
% virus = gseData.Data(:,contains(string(gseData.Header.Samples.source_name_ch1), 'VIRUS'));        % viral samples
% dataset_x = horzcat(bacteria, virus);
% dataset_y = horzcat(repmat(ones(1),1,size(bacteria, 2)) + 1,repmat(zeros(1),1,size(virus, 2)) + 1);
% % get most variable genes
% dsxp = dataset_x';
% dsxp_var = dsxp(:,var(dsxp) > 6000);
% dataset_x = dsxp_var';
% %%%

%%% LOAD BIOLOGY DATA
gseData = geoseriesread('/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE33341-GPL1261_series_matrix.txt');
sampleSources = unique(gseData.Header.Samples.source_name_ch1);
bacteria = gseData.Data(:,contains(string(gseData.Header.Samples.source_name_ch1), "S. aureus"));   % S.aureus samples
%bacteria = gseData.Data(:,(contains(string(gseData.Header.Samples.source_name_ch1), "S. aureus") & ~contains(string(gseData.Header.Samples.source_name_ch1), "2 hours")));
virus = gseData.Data(:,contains(string(gseData.Header.Samples.source_name_ch1), "E. coli"));        % E. colisamples
%virus = gseData.Data(:,(contains(string(gseData.Header.Samples.source_name_ch1), "E. coli") & ~contains(string(gseData.Header.Samples.source_name_ch1), "2 hours")));
dataset_x = horzcat(bacteria, virus);
dataset_y = horzcat(repmat(ones(1),1,size(bacteria, 2)) + 1,repmat(zeros(1),1,size(virus, 2)) + 1);

%playing around with the bootstrapping idea
%awgn(double(dataset_x(:,1)), 10)
%[bootstat_bacteria,bootsam_bacteria] = bootstrp(100,@mean,bacteria);
%[bootstat_virus,bootsam_virus] = bootstrp(100,@mean,virus);
%dataset_x = horzcat(bacteria, bootsam_bacteria, virus, bootsam_virus);
%dataset_y = horzcat(repmat(ones(1),1,size(bacteria, 2)) + 1, repmat(ones(1),1,size(bootsam_bacteria, 2)) + 1,repmat(zeros(1),1,size(virus, 2)) + 1,repmat(ones(1),1,size(bootsam_virus, 2)) + 1);


% get most variable genes
%dsxp = dataset_x';
%dsxp_var = dsxp(:,var(dsxp) > .5);
%size(dsxp_var)
%dataset_x = dsxp_var';
%%%


for i = 1:length(EXP_RANGE)
    
    flip_i = EXP_RANGE(i);
    
    for j = 1:length(EXP_RANGE_J)
        
        flip_j = EXP_RANGE_J(j);
        index = index + 1;
        
        "flip_i"
        flip_i
        
        "flip_j"
        flip_j
        
        for k = 1:ITER
            
            %[x, y, ff, xx, tt, dd] = genData(2,DIM,DS_SIZE,1000,1,'gen');
            [x, y, ff, xx, tt, dd] = subsetData(dataset_x, dataset_y, .8);

            Xt = double(x);
            yt = y;
            Xs = double(xx);
            ys = tt;
            
            % want to bootstrap new training samples here [bootstat,bootsam] = bootstrp(100,@mean,Xt);
            %bootsam_1 = bootstrap2(double( Xt(yt==1,:)'), 100);
            %bootsam_2 = bootstrap2(double( Xt(yt==2,:)'), 100);
            %
            %bootsam_1 = add_noise(double( Xt(yt==1,:)'), 100);
            %bootsam_2 = add_noise(double( Xt(yt==2,:)'), 100);        
            %Xt = horzcat(Xt', bootsam_1', bootsam_2')';
            %yt = horzcat(y', repmat(zeros(1), 1, size(bootsam_2,1)) + 1, repmat(ones(1), 1, size(bootsam_1,1)) + 1)';
            % FINISHED BOOTSTRAP SAMPLE
            
            [Xt, Xs] = standardise(Xt,Xs);  % get a "step direction is illegal" error if you don't standardize

            
            [yz, fdz] = injectLabelNoise(yt, [flip_i flip_j]);
            fprintf('[i = %f j = %f] sum=%f flip_ids=%f\n',flip_i,flip_j,sum(fdz ~= -1),(flip_i*length(yt) + flip_j*length(yt))/2);


            common_reg = 'lasso';
            common_sn = 1e-8;


            % add random/non-random label noise
            target = randi([1 CLS+1])-1;

            % winit and train a model prior to adding noise!
            %winit = randn(DIM+1,1);
            [idx2, ~] = rankfeatures(Xt', yt, 'Criterion', 'wilcoxon');  % features from known un-flipped, true data
            winit = randn(size(Xt(:,idx2(1:n_features)),2) + 1, 1);
            % rLR - standard log-reg w/o noise (for comparing results)
            options.estG = false;
            options.regFunc = common_reg;
            options.sn = common_sn;
            options.verbose = false;
            [w_nonoise, g_nonoise, l_lr_nonoise] = rlr(winit,eye(2),addbias(Xt(:,idx2(1:n_features))),yt,options);
            es_lr_nonoise(i,j,k) = sum(sign(addbias(Xs(:,idx2(1:n_features)))*w_nonoise) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs(:,idx2(1:n_features)))*w_nonoise,2);
            es_lr_nonoise_auc(i,j,k) = AUC;

            %size(Xt)
            sum(abs(w_nonoise) < 1e-6)/length(w_nonoise)

            
            % basic feature selection including all mislabelled samples
            if feature_select
                [idx, ~] = rankfeatures(Xt', yz, 'Criterion', 'wilcoxon');
                [idx2, ~] = rankfeatures(Xt', yt, 'Criterion', 'wilcoxon');  % features from known un-flipped, true data
                feature_lengths = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000];
                which_length = zeros(1,length(feature_lengths));
                for fl = 1:length(feature_lengths)
                    overlap = intersect(idx(1:feature_lengths(fl)), idx2(1:feature_lengths(fl)));
                    which_length(fl) = length(overlap) / feature_lengths(fl);
                end
                "which_length:"
                which_length
                Xt = Xt(:,idx(1:n_features)); % select top n_features features 
                Xs = Xs(:,idx(1:n_features)); % select top n_features features  
            end
            
            
            if use_PCs
                % run PCA on original dataset
                [COEFF, SCORE, explained_var] = pca(Xt);
                % select # of PCs to use in analysis
                n_pcs = sum((cumsum(explained_var)/sum(explained_var) < .9)); % use all PCs required to eplain 90% of variance                
                W = diag(std(Xs))\COEFF;  % project test dataset into PCA space
                [~, mu, we] = zscore(double(Xs));  % Getting mean and weights of data (for future data)
                we(we==0) = 1;
                xfs = double(Xs);  % New points in original feature space
                xfs = bsxfun(@minus, xfs, mu);
                xfs = bsxfun(@rdivide, xfs, we);
                projected_data = xfs*W; % New coordinates as principal components
                % set the data to PCs for both training and test sets
                Xt = SCORE(:,1:n_pcs);
                Xs = projected_data(:,1:n_pcs);
            end
            
            
            % create a new winit for training the noised, feature-selected !
            %winit = randn(DIM+1,1);
            winit = randn(size(Xt,2) + 1, 1);
            
           
            
            % rLR
            options.estG = false;
            options.regFunc = common_reg;
            options.sn = common_sn;
            options.verbose = false;
            [w, g, l_lr] = rlr(winit,eye(2),addbias(Xt),yz,options);
            es_lr(i,j,k) = sum(sign(addbias(Xs)*w) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*w,2);
            es_lr_auc(i,j,k) = AUC;
            
            sum(abs(w) < 1e-6)/length(w)
            
%             % rLR
%             options.estG = true;
%             options.regFunc = common_reg;
%             options.sn = common_sn;
%             options.verbose = false;
%             rr = .2;
%             [wr, gr, l_rlr] = rlr(winit,[1-rr rr;rr 1-rr],addbias(Xt),yz,options);
%             es_rlr(i,j,k) = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
%             [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr,2);
%             es_rlr_auc(i,j,k) = AUC;
            
            % rLR
            options.estG = true;
            options.regFunc = common_reg;
            options.sn = common_sn;
            options.verbose = false;
            rr = .2;
            [wr, gr, l_rlr] = rlr(winit,[1-rr rr;rr 1-rr],addbias(Xt),yz,options);
            es_rlr(i,j,k) = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr,2);
            es_rlr_auc(i,j,k) = AUC;
           
            sum(abs(wr) < 1e-6)/length(wr)
            
            %rLR iteration 2 after filtering genes
            keep = abs(wr) > 1e-6;
            keep = keep(1:length(keep)-1);
            options.estG = true;
            options.regFunc = common_reg;
            options.sn = common_sn;
            options.verbose = false;
            winit2 = randn(sum(keep) + 1, 1);
            rr = .2;
            [wr2, gr2, l_rlr2] = rlr(winit2,[1-rr rr;rr 1-rr],addbias(Xt(:,keep)),yz,options);
            es_rlr2(i,j,k) = sum(sign(addbias(Xs(:,keep))*wr2) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs(:,keep))*wr2,2);
            es_rlr_auc2(i,j,k) = AUC;
            sum(abs(wr2) < 1e-6)/length(wr2)
            
            % nLR
            options.estG = true;
            options.regFunc = common_reg;
            options.sn = common_sn;
            options.verbose = false;
            [wg, nd, l_gam] = gammalr(winit,addbias(Xt), yz, options);
            es_gammalr(i,j,k) =  sum(sign(addbias(Xs)*wg) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wg,2);
            es_gammalr_auc(i,j,k) = AUC;
            
            sum(abs(wg) < 1e-6)/length(wg)
            
            %fprintf('AUC lr=%f rlr1=%f rlr2=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_rlr_auc2(i,j,k),es_gammalr_auc(i,j,k));
            fprintf('AUC lr=%f rlr1=%f rlr2=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_rlr_auc2(i,j,k),es_gammalr_auc(i,j,k));
            fprintf('ERR lr=%f rlr1=%f rlr2=%f gammalr=%f\n',es_lr(i,j,k),es_rlr(i,j,k),es_rlr2(i,j,k),es_gammalr(i,j,k));
            
            
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
plot_heatmap(es_rlr2, 'MEAN ES\_RLR_noGest', 'figure2b', 0, .55, EXP_RANGE_J, EXP_RANGE)
%plot_heatmap(es_rlr3, 'MEAN ES\_RLR\_Gest', 'figure2c', 0, .55, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr, 'MEAN ES\_GAMMALR', 'figure3', 0, .55, EXP_RANGE_J, EXP_RANGE)

plot_heatmap(es_lr_auc, 'MEAN ES\_LR AUC', 'figure1_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_rlr_auc, 'MEAN ES\_RLR AUC', 'figure2_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_rlr_auc2, 'MEAN ES\_RLR\_noGest AUC', 'figure2b_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
%plot_heatmap(es_rlr_auc3, 'MEAN ES\_RLR\_Gest AUC', 'figure2c_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr_auc, 'MEAN ES\_GAMMALR AUC', 'figure3_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)

plot_lineplot(es_lr, es_rlr, es_gammalr ,EXP_RANGE_J, EXP_RANGE, 'Figure4')
plot_lineplot(100 - es_lr, 100 - es_rlr, 100 - es_gammalr ,EXP_RANGE_J, EXP_RANGE, 'Figure5')
plot_lineplot(es_lr_auc, es_rlr_auc, es_gammalr_auc ,EXP_RANGE_J, EXP_RANGE, 'Figure5_auc')












colormap(winter)
scatter(Xt(:,1), Xt(:,2), 100, yt, 's', 'filled');
hold on
scatter(Xs(:,1), Xs(:,2), 100, ys, 'o', 'filled');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure7','-bestfit','-dpdf')
hold off;

colormap(winter)
scatter(Xt(:,1), Xt(:,2), 100, yz, 's', 'filled');
hold on
scatter(Xs(:,1), Xs(:,2), 100, ys, 'o', 'filled');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure7b','-bestfit','-dpdf')
hold off;


% PCA projection of training data (original PCA)
[COEFF,SCORE] = pca(Xt);
scatter(SCORE(:,1), SCORE(:,2));
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure6','-bestfit','-dpdf')


% PCA projection of test data
W = diag(std(Xs))\COEFF;
% Getting mean and weights of data (for future data)
[~, mu, we] = zscore(Xs);
we(we==0) = 1;
% New points in original feature space
xfs = Xs;
xfs = bsxfun(@minus, xfs, mu);
xfs = bsxfun(@rdivide, xfs, we);
% New coordinates as principal components
projected_data = xfs*W;
% plot the new projection on top of original PCA space
scatter(SCORE(:,1), SCORE(:,2));
hold on;
scatter(y(:,1), y(:,2), 'red');






(es_rlr_auc - es_lr_auc) > 0
%(es_rlr_auc2 - es_lr_auc) > 0
%(es_rlr_auc3 - es_lr_auc) > 0