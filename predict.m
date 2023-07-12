function [ztrain,ztest,zaug] = predict(train_data,test_data,nameMethod,lam)
%trains model on train data
%train_data = training data
%test_data = test data
%nameMethod = name of the denoising method to use
%ztrain = denoised training data
%ztest = interpolated test data

% SASS
if strcmp(nameMethod,'SASS')
    sorttrainData = sortrows(train_data,1);
    y = sorttrainData(:,2);
    
    d = 2;          % d : filter order parameter (d = 1, 2, or 3)
    fc = 0.05;     % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
    K = 2;          % K : order of sparse derivative
    %lam = 1.2;      % lam : regularization parameter
    
    [ztrain, ~] = sass_L1(y, d, fc, K, lam);
    
    ztest = interp1(sorttrainData(:,1),ztrain,test_data(:,1),"pchip",'extrap');
    
    augdata = [[sort(train_data(:,1)) ztrain];[test_data(:,1) ztest]]; 
    sortedaugdata = sortrows(augdata,1);
    zaug = sortedaugdata(:,2); 
end
%BEADS
if strcmp(nameMethod,'BEADS')
    sorttrainData = sortrows(train_data,1);
    y = sorttrainData(:,2);
    
    fc = 0.006;             % fc : cut-off frequency (cycles/sample)
    d = 1;                  % d : filter order parameter (d = small positive integer)
    Nit = 50;
    EPS0 = 1e-5;            % for x itself
    EPS1 = 1e-5;            % for derivatives
    %alpha = 1;
    
    r = 1;                  % to set the asymmetric ratio
    pen = 'L1_v1';          % to choose penalty function
    
    [x,f,~,~,~] = beads(y, d, fc, r, lam, pen, Nit, EPS0, EPS1);
    
    ztrain = x+f; %taking accound of baseline
    ztest = interp1(sorttrainData(:,1),ztrain,test_data(:,1),"pchip",'extrap');
    
    augdata = [[sort(train_data(:,1)) ztrain];[test_data(:,1) ztest]]; 
    sortedaugdata = sortrows(augdata,1);
    zaug = sortedaugdata(:,2); 
end
end