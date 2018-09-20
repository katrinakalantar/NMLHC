% %==========================================================================
% % Name  : A function for evaluating a binary logistic regression classifier
% % Author: Jakramate Bootkrajang
% % Input : x  = input data
% %         y  = target values
% %         w  = learned parameters vector of the LR
% %         g  = the label flipping matrix
% %         pyHat = predicting the observed label ? 
% % Return: pr  = probability of predicted class
% %         pd  = predicted class  {-1,1}
% %         e   = misclassification vector
% %         er  = misclassification rate
% %==========================================================================
% %==========================================================================
% 
% function [pr pd e er] = evalLR(x, y, w, g, pyHat)
% 
% % the function uses {-1,1} class label
% %y = castLabel(y,-1);
% 
% t = (x*w);
% if pyHat   
%     % calculating posterior probability of the 'observed' positive samples
%     pPos = (g(1,2) * (1 ./ ((1 ./ exp(-t))+1))) + (g(2,2) ./ (1+exp(-t)));
% else
%     % calculating posterior probability of the 'true' positive samples 
%     pPos = sigmoid(t);         
% end
% pPos    = [1-pPos pPos];
% [pr pd] = max(pPos, [], 2);
% %pd      = castLabel(pd,-1);  
% e       = (pd ~= y);
% er      = sum(e)/size(x,1);



% Function to evaluate model of logistic regression family
% Author: Jakramate Bootkrajang
% Input : x  = Samples to be classified
%         y  = Target value
%         w  = Weight vector for Sigmoid/Softmax
%         fd = Label flipping indicator
% Return: prob  = probability of predicted class
%         pred  = predicted class
%         eIdx  = misclassification vector
%         eRate = misclassification rate


function [prob pred eIdx eRate] = evalLR(x, y, w)

if (size(x,2) ~= size(w,1))
    error('input and weight vector have different length missing bias?');
end

TPNT = size(x,1);                % total number of test samples 
CLS  = size(w,2) + 1;            % determine number of classes

if (CLS == 2)        
    postY  = sigmoid(x * w);     % calculate the sigmoid
    postY  = [1-postY postY]';  
    %postYh = gamma' * postY;    % used to predict given labels
%     tmp    = gamma(:,yz) .* postY;
%     Pn     = tmp(2,:) ./ sum(tmp,1);        
else   
    postY = softmax(x * w)';     % calculate softmax for prediction        
end

[prob pred] = max(postY); 
eIdx        = (pred' ~= y);
eRate       = sum(eIdx)/TPNT;

% pred2       = (Pn > 0.5) + 1;
% eIdx2       = (pred2' ~= y);
% eRate2      = (sum(eIdx2)/TPNT) * 100;
