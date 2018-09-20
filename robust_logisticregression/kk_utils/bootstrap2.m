function bootstrap_sample = bootstrap2(dataset, n_boot)
	fprintf("hlelo")
    R = dataset';
    bootstrap_sample = zeros(n_boot, size(R,2));
    for i = 1:n_boot
        samples = R(sparse(randi(size(R,1),size(R,2),1), 1:size(R,2), true));
        bootstrap_sample(i,:) = samples';
    end
end




