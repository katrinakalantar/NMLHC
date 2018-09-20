function cK = compositeKernel(beta, k)

cK = zeros(size(k{1}));
p  = length(beta);

for i=1:p
    cK = cK + beta(i) * k{i};
end

cK(:,1) = 1;  % adjust bias
