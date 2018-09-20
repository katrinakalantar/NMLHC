function noised_sample = add_noise(dataset, n_boot)
    R = dataset';
    noised_sample = zeros(n_boot, size(R,2));
    for i = 1:n_boot
        current_sample = randsample([1:size(R,1)],1);
        noised_sample(i,:) = awgn(R(current_sample,:), 10);
    end
end