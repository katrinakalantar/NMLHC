function [dataset_train, dataset_test] = override_kmeans(dataset_train, dataset_test, k, yval,ytrue, run_pca)

    tic
    
    %goodsamples = zeros(1,size(dataset_train,1)); % create vector of zeros for each sample

    if run_pca
        [~, SCORE, explained_var] = pca(dataset_train);
        n_pcs = sum((cumsum(explained_var)/sum(explained_var) < .9));
        E = evalclusters(SCORE(:,1:n_pcs),'kmeans','silhouette','klist',[1:6]);
        E
        k_means_result = kmeans(SCORE(:,1:n_pcs), E.OptimalK);
    else
        E = evalclusters(dataset_train,'kmeans','silhouette','klist',[1:6]);
        E
        k_means_result = kmeans(dataset_train,E.OptimalK); % get k-means result
    end
    
     total_ones = sum(yval == 1);
     total_twos = sum(yval == 2);
    
    
    mat1 = [k_means_result, yval];
    
    
    cluster_ids = zeros(1, length(E.InspectedK));
    good_samples = zeros(1, length(yval));
    
    for curr_kval = 1:E.OptimalK
        [a1,b1] = hist(mat1(mat1(:,1)==curr_kval,2), unique(mat1(mat1(:,1)==curr_kval,2)));
        a1 = a1 ./ [total_ones, total_twos];
        if max(a1)/sum(a1) > .5   
            cluster_ids(curr_kval) =  b1(a1 == max(a1));
        else
            cluster_ids(curr_kval) = 0;
        end        
        fprintf('--> k = %f;  A_1 = %f, A_2 = %f ;   percent = %f, clusterID = %f\n', ...
            curr_kval, a1(1), a1(2), max(a1)/sum(a1), cluster_ids(curr_kval))
        
        samples_in_cluster = (k_means_result  == cluster_ids(curr_kval)); 
        samples_correctly_in_cluster = (yval == cluster_ids(curr_kval)) &  samples_in_cluster;
        
        good_samples = good_samples | samples_correctly_in_cluster';
        
        fprintf('--> samples_correctly_in_cluster= %f  \n', sum(samples_correctly_in_cluster)/sum(samples_in_cluster))
        
    end

    ds_train = dataset_train(good_samples==1,:);
    yval_train = yval(good_samples == 1);
   
    
    % Use rankfeatures function to get significant features all at once
    [idx, ~] = rankfeatures(ds_train', yval_train, 'Criterion', 'wilcoxon');
    
    %%% this is all for benchmarking purposes (which features are selected?)
    [idx_full, ~] = rankfeatures(dataset_train', yval, 'Criterion', 'wilcoxon');
    [idx_full2, ~] = rankfeatures(dataset_train', ytrue, 'Criterion', 'wilcoxon');  % features from known un-flipped, true data
    length(idx)
    length(idx_full)
    length(idx_full2)
    i1 = length(intersect(idx(1:5000), idx_full(1:5000)))/5000;
    i2 = length(intersect(idx(1:5000), idx_full2(1:5000)))/5000;
    i3 = length(intersect(idx_full(1:5000), idx_full2(1:5000)))/5000;
    fprintf('--> INTERSECT fancy v. basic mislabelled:%f\n', i1)
    fprintf('--> INTERSECT fancy v. optimal:%f\n', i2);
    fprintf('--> INTERSECT basic mislabelled v. optimal:%f\n', i3);
    %%% end benchmarking code
    
    
    dataset_train = dataset_train(:,idx(1:1000)); % select top 1000 features 
    dataset_test = dataset_test(:,idx(1:1000)); % select top 1000 features 
    
    
%     % Use anova for selecting, iterate over each column; compute anova; get
%     % the significant features; TAKES FOREVER TO RUN.
%     tt = zeros(1, size(ds_train,2));
%     for ind = 1:size(ds_train,2)
%         %tt(i) = ttest2(Xt(:,i), yt);
%         [p, ~, ~] = anova1(ds_train(:,ind),yval_train, 'off');  % yz' or yt
%         tt(ind) = p;
%     end
%     sxt = size(ds_train,2);
%     dataset_train = dataset_train(:,(tt < .0001/sxt));
% 	dataset_test = dataset_test(:,(tt < .0001/sxt));

    toc

end