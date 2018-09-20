function x = addbias(x)

nrow = size(x,1);
x    = [ones(nrow,1) x];