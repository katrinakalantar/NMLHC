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
ITER            = 5; %10;%10; 
options.maxIter = 100;

% parameters iterated
%EXP_RANGE       = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5]; %y-axis of heatmap
%EXP_RANGE_J     = [0, .01, .025, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap
EXP_RANGE       = [0, .05, .1, .15, .2, .25, .3];%, .05, .1, .15, .2, .3, .4, .5]; %y-axis of heatmap.025, 
EXP_RANGE_J     = [0, .05, .1, .15, .2, .25, .3];%, .05, .1, .15, .2, .3, .4, .5];  % x-axis of heatmap

CLS             = 2;    % number of classes
DIM             = 1000; % dimensionality of dataset generated
DS_SIZE         = 50;  % dataset size

% preallocating error storage
es_lr      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_lr_nonoise      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
% es_rlr2    = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
% es_rlr3    = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_gammalr = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);

% preallocating AUC storage
es_lr_auc       = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
es_rlr_auc      = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
% es_rlr_auc2     = nan(length(EXP_RANGE), length(EXP_RANGE_J), ITER);
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

%             tic
%             if feature_select         
%                 % trying out t-test for feature selection
%                 tt = zeros(1, size(Xt,2));
%                 for ind = 1:size(Xt,2)
%                     %tt(i) = ttest2(Xt(:,i), yt);
%                     [p, t, stats] = anova1(Xt(:,ind),yz', 'off');  % yz' or yt
%                     tt(ind) = p;
%                 end
%                 sxt = size(Xt,2);
%                 Xt = Xt(:,(tt < .0001/sxt));
%                 Xs = Xs(:,(tt < .0001/sxt));
%             end
%             toc

            
            %size(Xt)
        
            
            % train a gold standard decision boundary for noise generation
            % using LR
            options.estG    = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            %[ww, ~, ~]      = rlr(winit,eye(CLS),addbias(Xt),yt,options);
            
            [yz, fdz] = injectLabelNoise(yt, [flip_i flip_j]);
            fprintf('[i = %f j = %f] sum=%f flip_ids=%f\n',flip_i,flip_j,sum(fdz ~= -1),(flip_i*length(yt) + flip_j*length(yt))/2);
            
%             % fancy function for deriving features from only "robust"
%             % samples
%             if feature_select
%                [Xt, Xs] = override_kmeans(Xt, Xs, 2, yz, yt, true);
%             end
            
            % size(Xt)
            
%             tic
%             if feature_select         
%                 % trying out t-test for feature selection
%                 tt = zeros(1, size(Xt,2));
%                 for ind = 1:size(Xt,2)
%                     %tt(i) = ttest2(Xt(:,i), yt);
%                     [p, t, stats] = anova1(Xt(:,ind),yt, 'off'); %yz'
%                     tt(ind) = p;
%                 end
%                 sxt = size(Xt,2);
%                 Xt = Xt(:,(tt < .0001/sxt));
%                 Xs = Xs(:,(tt < .0001/sxt));
%             end
%             toc



            % add random/non-random label noise
            target = randi([1 CLS+1])-1;

            % winit and train a model prior to adding noise!
            %winit = randn(DIM+1,1);
            [idx2, ~] = rankfeatures(Xt', yt, 'Criterion', 'wilcoxon');  % features from known un-flipped, true data
            winit = randn(size(Xt(:,idx2(1:n_features)),2) + 1, 1);
            % rLR - standard logistic regression with no noise (for
            % comparing results)
            options.estG = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            [w_nonoise, g_nonoise, l_lr_nonoise] = rlr(winit,eye(2),addbias(Xt(:,idx2(1:n_features))),yt,options);
            es_lr_nonoise(i,j,k) = sum(sign(addbias(Xs(:,idx2(1:n_features)))*w_nonoise) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs(:,idx2(1:n_features)))*w_nonoise,2);
            es_lr_nonoise_auc(i,j,k) = AUC;


            
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
                which_length
                Xt = Xt(:,idx(1:n_features)); % select top n_features features 
                Xs = Xs(:,idx(1:n_features)); % select top n_features features  
            end
            
            
            if use_PCs
                % run PCA on original dataset
                [COEFF, SCORE, explained_var] = pca(Xt);
                % select # of PCs to use in analysis
                n_pcs = sum((cumsum(explained_var)/sum(explained_var) < .9)); % use all PCs required to eplain 70% of variance                
                % project test dataset into PCA space
                W = diag(std(Xs))\COEFF;
                [~, mu, we] = zscore(double(Xs));  % Getting mean and weights of data (for future data)
                we(we==0) = 1;
                xfs = double(Xs);  % New points in original feature space
                xfs = bsxfun(@minus, xfs, mu);
                xfs = bsxfun(@rdivide, xfs, we);
                projected_data = xfs*W; % New coordinates as principal components
                % set the data to PC representation for both training and
                % test sets
                Xt = SCORE(:,1:n_pcs);
                Xs = projected_data(:,1:n_pcs);
            end
            
            
            % create a new winit for training the noised, feature-selected !
            %winit = randn(DIM+1,1);
            winit = randn(size(Xt,2) + 1, 1);
            
            
            % rLR
            options.estG = false;
            options.regFunc = 'lasso';
            options.verbose = false;
            [w, g, l_lr] = rlr(winit,eye(2),addbias(Xt),yz,options);
            es_lr(i,j,k) = sum(sign(addbias(Xs)*w) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*w,2);
            es_lr_auc(i,j,k) = AUC;
            
            % rLR
            options.estG = true;
            options.regFunc = 'lasso';
            options.verbose = false;
            rr = .2;
            [wr, gr, l_rlr] = rlr(winit,[1-rr rr;rr 1-rr],addbias(Xt),yz,options);
            es_rlr(i,j,k) = sum(sign(addbias(Xs)*wr) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr,2);
            es_rlr_auc(i,j,k) = AUC;
            
%             % rLR without estimating G, but with better estimation of G up
%             % front
%             options.estG = false;
%             options.regFunc = 'lasso';
%             options.verbose = false;
%             [wr2, gr2, l_rlr2] = rlr(winit,[1-flip_i flip_i;flip_j 1-flip_j], addbias(Xt), yz, options);
%             es_rlr2(i,j,k) = sum(sign(addbias(Xs)*wr2) ~= castLabel(ys,-1))/length(ys);
%             [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr2,2);
%             es_rlr_auc2(i,j,k) = AUC;
%             
%             % rLR with estimating G AND better estimation of G up front
%             options.estG = true;
%             options.regFunc = 'lasso';
%             options.verbose = false;
%             [wr3, gr3, l_rlr3] = rlr(winit,[1-flip_i flip_i;flip_j 1-flip_j], addbias(Xt), yz, options);
%             es_rlr3(i,j,k) = sum(sign(addbias(Xs)*wr3) ~= castLabel(ys,-1))/length(ys);
%             [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wr3,2);
%             es_rlr_auc3(i,j,k) = AUC;
           
            % nLR
            options.estG = true;
            options.regFunc = 'lasso';
            options.verbose = false;
            [wg, nd, l_gam] = gammalr(winit,addbias(Xt), yz, options);
            es_gammalr(i,j,k) =  sum(sign(addbias(Xs)*wg) ~= castLabel(ys,-1))/length(ys);
            [~,~,~,AUC] = perfcurve(ys,addbias(Xs)*wg,2);
            es_gammalr_auc(i,j,k) = AUC;
            
            %fprintf('AUC lr=%f rlr1=%f rlr2=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_rlr_auc2(i,j,k),es_gammalr_auc(i,j,k));
            fprintf('AUC lr=%f rlr1=%f gammalr=%f\n',es_lr_auc(i,j,k),es_rlr_auc(i,j,k),es_gammalr_auc(i,j,k));
            
            
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
%plot_heatmap(es_rlr2, 'MEAN ES\_RLR_noGest', 'figure2b', 0, .55, EXP_RANGE_J, EXP_RANGE)
%plot_heatmap(es_rlr3, 'MEAN ES\_RLR\_Gest', 'figure2c', 0, .55, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr, 'MEAN ES\_GAMMALR', 'figure3', 0, .55, EXP_RANGE_J, EXP_RANGE)

plot_heatmap(es_lr_auc, 'MEAN ES\_LR AUC', 'figure1_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_rlr_auc, 'MEAN ES\_RLR AUC', 'figure2_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
%plot_heatmap(es_rlr_auc2, 'MEAN ES\_RLR\_noGest AUC', 'figure2b_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
%plot_heatmap(es_rlr_auc3, 'MEAN ES\_RLR\_Gest AUC', 'figure2c_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)
plot_heatmap(es_gammalr_auc, 'MEAN ES\_GAMMALR AUC', 'figure3_auc', 0, 1, EXP_RANGE_J, EXP_RANGE)




%[COEFF,SCORE] = pca(Xt)


x = EXP_RANGE_J;
x2 = EXP_RANGE;

y = mean(es_lr,3).*100;
err = std(es_lr,0,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
plot(x,y(1,:), 'LineWidth', 2,'color',[89 47 188]./255);  hold on;
e = errorbar(x,y(1,:),err(1,:),'-k.');
e.Color = [89 47 188]./255; %purple
hold on;
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
plot(x2,y(:,1), 'LineWidth', 2,'color', [165 127 255]./255);  hold on;
e = errorbar(x2,y(:,1),err(:,1));
e.Color = [165 127 255]./255; %light purple
hold on;

y = mean(es_rlr,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_rlr,0,3).*100;
plot(x,y(1,:), 'LineWidth', 2,'color',[40 155 71]./255);  hold on;
e2 = errorbar(x,y(1,:),err(1,:));
e2.Color = [40 155 71]./255; %dark green
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_rlr,0,3).*100;
plot(x2,y(:,1), 'LineWidth', 2,'color', [53 214 96]./255);  hold on;
e2 = errorbar(x2,y(:,1),err(:,1));
e2.Color = [53 214 96]./255; %light green
hold on;

y = mean(es_gammalr,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_gammalr,0,3).*100;
plot(x,y(1,:), 'LineWidth', 2,'color',[47 170 188]./255);  hold on;
e3 = errorbar(x,y(1,:),err(1,:));
e3.Color = [47 170 188]./255; % dark blue
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_gammalr,0,3).*100;
plot(x2,y(:,1), 'LineWidth', 2,'color', [60 218 242]./255);  hold on;
e3 = errorbar(x2,y(:,1),err(:,1));
e3.Color = [60 218 242]./255; % light blue
hold on;

xlabel("% Flipped")
ylabel("Error")
ylim([0 100])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure4','-bestfit','-dpdf')

hold off;


x = EXP_RANGE_J;
x2 = EXP_RANGE;

y = 100 - mean(es_lr,3).*100;
err = std(es_lr,0,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
plot(x,y(1,:), 'LineWidth', 2,'color',[89 47 188]./255);  hold on;
e = errorbar(x,y(1,:),err(1,:),'-k.');
e.Color = [89 47 188]./255; %purple
hold on;
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
plot(x2,y(:,1), 'LineWidth', 2,'color', [165 127 255]./255);  hold on;
e = errorbar(x2,y(:,1),err(:,1));
e.Color = [165 127 255]./255; %light purple
hold on;

y = 100 - mean(es_rlr,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_rlr,0,3).*100;
plot(x,y(1,:), 'LineWidth', 2,'color',[40 155 71]./255);  hold on;
e2 = errorbar(x,y(1,:),err(1,:));
e2.Color = [40 155 71]./255; %dark green
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_rlr,0,3).*100;
plot(x2,y(:,1), 'LineWidth', 2,'color', [53 214 96]./255);  hold on;
e2 = errorbar(x2,y(:,1),err(:,1));
e2.Color = [53 214 96]./255; %light green
hold on;

y = 100 - mean(es_gammalr,3).*100;
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_gammalr,0,3).*100;
plot(x,y(1,:), 'LineWidth', 2,'color',[47 170 188]./255);  hold on;
e3 = errorbar(x,y(1,:),err(1,:));
e3.Color = [47 170 188]./255; % dark blue
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_gammalr,0,3).*100;
plot(x2,y(:,1), 'LineWidth', 2,'color', [60 218 242]./255);  hold on;
e3 = errorbar(x2,y(:,1),err(:,1));
e3.Color = [60 218 242]./255; % light blue
hold on;


xlabel("% Flipped")
ylabel("Accuracy (%)")
ylim([0 100])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure5','-bestfit','-dpdf')
hold off;



x = EXP_RANGE_J;
x2 = EXP_RANGE;

y = mean(es_lr_auc,3);
err = std(es_lr_auc,0,3);
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
plot(x,y(1,:), 'LineWidth', 2,'color',[89 47 188]./255);  hold on;
e = errorbar(x,y(1,:),err(1,:),'-k.');
e.Color = [89 47 188]./255; %purple
hold on;
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
plot(x2,y(:,1), 'LineWidth', 2,'color', [165 127 255]./255);  hold on;
e = errorbar(x2,y(:,1),err(:,1));
e.Color = [165 127 255]./255; %light purple
hold on;

y = mean(es_rlr_auc,3);
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_rlr_auc,0,3);
plot(x,y(1,:), 'LineWidth', 2,'color',[40 155 71]./255);  hold on;
e2 = errorbar(x,y(1,:),err(1,:));
e2.Color = [40 155 71]./255; %dark green
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_rlr_auc,0,3);
plot(x2,y(:,1), 'LineWidth', 2,'color', [53 214 96]./255);  hold on;
e2 = errorbar(x2,y(:,1),err(:,1));
e2.Color = [53 214 96]./255; %light green
hold on;

y = mean(es_gammalr_auc,3);
% mean/error along x-axis of heatmap // % Flipped in J, one to zero
err = std(es_gammalr_auc,0,3);
plot(x,y(1,:), 'LineWidth', 2,'color',[47 170 188]./255);  hold on;
e3 = errorbar(x,y(1,:),err(1,:));
e3.Color = [47 170 188]./255; % dark blue
% mean/error along y-axis of heatmap // % Flipped in I, zero to one
err = std(es_gammalr_auc,0,3);
plot(x2,y(:,1), 'LineWidth', 2,'color', [60 218 242]./255);  hold on;
e3 = errorbar(x2,y(:,1),err(:,1));
e3.Color = [60 218 242]./255; % light blue
hold on;


xlabel("% Flipped")
ylabel("AUC")
ylim([0 1])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('figure5_auc','-bestfit','-dpdf')
hold off;









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






%TESTING

% create random matrix
a = 0; b = 1:10;
B = repmat(b,5,1);
R = unifrnd(a,B); % 5 samples (rows), 10 features (cols).

R2 = bootstrap2(R', 4);  % it does look like each sampled value comes from the appropriate column (features)

size(Xt) %130 samples (rows), 3000 features (cols)
bootsam_1 = bootstrap2(double( Xt(yt==1,:)'), 100);  % 
bootsam_2 = bootstrap2(double( Xt(yt==2,:)'), 100);
            


(es_rlr_auc - es_lr_auc) > 0
%(es_rlr_auc2 - es_lr_auc) > 0
%(es_rlr_auc3 - es_lr_auc) > 0