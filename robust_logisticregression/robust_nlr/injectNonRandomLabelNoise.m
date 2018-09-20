% Name  : A function to inject non random label noise
% Author: Jakramate Bootkrajang
% Input : y           = true label {1,...n} representation
%         targetClass = {0,1,2....} s.t. 0 is symmetric NoiseRateing
%         NoiseRate   = (0,1)
% Output: yz          = noisy label  {1,...n} representation
%         fp          = flip indicator vector


function [yz fd] = injectNonRandomLabelNoise(y, fd, x, w, target, noisepdf, p0, p1)

yz  =  castLabel(y,-1);
y   =  castLabel(y,2);

% distance given w
dis = (x*w)/norm(w);

Z   =  abs(noisepdf(dis, p0, p1));

% sampling some numbers
if target == 0
    idx = find(rand(size(y)) <= Z);
    yz(idx) = yz(idx) * -1;
    fd(idx) = fd(idx) * -1;
else
    idx = find(rand(size(y)) <= Z & (y==target));
    yz(idx) = yz(idx) * -1;
    fd(idx) = fd(idx) * -1;
end

yz = castLabel(yz,2);