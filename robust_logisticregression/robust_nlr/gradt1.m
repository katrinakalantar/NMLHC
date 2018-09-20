function [fv dfv] = gradt1(t1, z, t, y, t0, k0, k1, options)

g01 = gammapdf(z,k0,t0);
g10 = gammapdf(z,k1,t1);
% g01 = gammapdf(z,k0,t0) .* (z~=0);
% g10 = gammapdf(z,k1,t1) .* (z~=0);
g00 = 1 - g01;
g11 = 1 - g10;

p1 = 1./ (1+exp(-t));
p0 = 1-p1;

% computing function value
s0    = (g00 .* p0) + (g10 .* p1);
s1    = (g01 .* p0) + (g11 .* p1);

s0(s0==0) = eps; 
s1(s1==0) = eps;

fv = -sum((y .* log(s1)) + ((1 - y).* log(s0)),1)...
     + options.regV + options.bReg0 * log(t0-1) + options.bReg1*log(t1-1);

% computing corresponding derivative
tmp1  = ((1 - y) ./ s0) - (y ./ s1);


tmp2  = ((g10.*z./t1^2) - (g10 ./t1));

gAux1 = (tmp1 .* tmp2) .* p1;

dfv   = -sum(gAux1,1)' + options.bReg1/log(t1-1);