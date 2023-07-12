function [z,b] = CVBEADS(y,alphas,K,type)
% denoising and baseline correction using BEADS with method of Cross-Validation for alpha tuning 
% Input
% y : noisy signal
% alphas : set of relation parameters
% K : number of cross-validation folds
% Output 
% z : denoised signal
% b : baseline

x=linspace(1,length(y),length(y))';
dataframe=[x y];

% BEADS fixed hyperparameters 
fc = 0.006;             % fc : cut-off frequency (cycles/sample)
d = 1;                  % d : filter order parameter (d = small positive integer)
Nit = 50;
EPS0 = 1e-5;            % for x itself
EPS1 = 1e-5;            % for derivatives
    
r = 1;                  % to set the asymmetric ratio
pen = 'L1_v1';          % to choose penalty function

% Cross-validation
conf_scores = zeros([2, length(alphas)]);
for i = 1:length(alphas)
alpha = alphas(i);
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
    
[ztrain,btrain,~,~,~] = beads(ytrain, d, fc, r, alpha, pen, Nit, EPS0, EPS1);
ztrain = ztrain + btrain;
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
opt_alpha = alphas(ind);

% compute optimal signal
[z,b,~,~,~] = beads(y, d, fc, r, opt_alpha, pen, Nit, EPS0, EPS1);

end