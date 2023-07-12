function [train_data,test_data] = splitdata(folds,fold)
%splits data into train and test data 
%given fold is test data, remaining train data
%folds=matrix of x,y,k where k index of fold
%fold=integer between 1 and k where k is number of folds
indTest = folds(:,3) == fold;
indTrain = folds(:,3) ~= fold;
test_data = folds(indTest,1:2);
train_data = folds(indTrain,1:2);
end