function [z] = CVSG(y,ws, K, type)
% denoising using SG with method of Cross-Validation for window size tuning 
% Input
% y : noisy signal
% lambdas : set of regularization parameters
% K : number of cross-validation folds
% Output 
% z : denoised signal

x=linspace(1,length(y),length(y))';
dataframe=[x y];

% SG fixed hyperparameters 
order = 2;

% Cross-validation
conf_scores = zeros([2, length(ws)]);
for i = 1:length(ws)
w = ws(i);
% shuffle
shuffledData = shuffle(dataframe);
% split
folds = split(shuffledData,K);
% split data into train and test set
scores = zeros([1, K]);
for k = 1:K
[train_data,test_data]=splitdata(folds,k);

% predict signal value for test set
sorttrainData = sortrows(train_data,1);
ytrain = sorttrainData(:,2);
    
ztrain = sgolayfilt(ytrain,order,w);
    
ztest = interp1(sorttrainData(:,1),ztrain,test_data(:,1),"pchip",'extrap');
    
augdata = [[sort(train_data(:,1)) ztrain];[test_data(:,1) ztest]]; 
sortedaugdata = sortrows(augdata,1);
zaug = sortedaugdata(:,2); 

% evaluate method on test data
score = evaluate(ztest,test_data(:,2),type);
scores(k) = score;
end
scores_mean = mean(scores);
scores_std = std(scores);
conf_scores(1,i) = scores_mean;
conf_scores(2,i) = scores_std;
end

% find best lambda
[min_CVerr,ind] = min(conf_scores(1,:));
opt_w = ws(ind); % lam : regularization parameter

% compute optimal signal
z = sgolayfilt(y,order,opt_w);

end