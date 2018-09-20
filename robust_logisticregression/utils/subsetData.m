% Author: Jakramate Bootkrajang
% Desc: Generate synthetic data generatively/discriminatively
% label from this function is 1...N

function [X Y FD Xt Yt FDt] = subsetData(INPUT_DATA, INPUT_LABELS, PCT)

ds = datasample([1:size(INPUT_DATA,2)], size(INPUT_DATA,2), 'Replace',false);
training = ds(1:round(length(ds)*PCT));
testing = ds(round(length(ds)*.8) + 1:length(ds));

X = INPUT_DATA(:,training)';
Y = INPUT_LABELS(training)';
FD = zeros(length(training),1);

Xt = INPUT_DATA(:,testing)';
Yt = INPUT_LABELS(testing)';
FDt = zeros(length(testing),1);

end

