%==========================================================================
% Name  : A function for evaluating a multinomial logistic regression classifier
% Author: Jakramate Bootkrajang
% Input : x  = input data
%         y  = target values
%         w  = learned parameters vector of the LR
%         g  = the label flipping matrix
%         pyHat = predicting the observed label ? 
% Return: pr  = probability of predicted class
%         pd  = predicted class  {-1,1}
%         e   = misclassification vector
%         er  = misclassification rate
%==========================================================================
%==========================================================================

function [pr pd e er] = evalmLR(x, y, w, g, pyHat)

% the function uses {-1,1} class label
%y = castLabel(y,2);

t = (x*w);
if pyHat   
    % calculating posterior probability of the 'observed' positive samples
    %pPos = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));
else
    % calculating posterior probability of the 'true' positive samples 
    pPos = softmax(t);         
end
[pr pd] = max(pPos, [], 2);  
e       = (pd ~= y);
er      = sum(e)/size(x,1);