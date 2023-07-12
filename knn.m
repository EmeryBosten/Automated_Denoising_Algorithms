function [ns] = knn(X,q,k)
% k-nearest neighbors
% Input
% X = dataset consisting of M datapoints
% q = query point
% k = positive integer representing number of nearest neighbors to find
% Output
% ns = a list of k nearest neighbors to q in X

[row,col] = size(X);
% Compute distance between q and each point in X
dist = zeros([row, 1]);
for i = 1:row
    dist(i) = euclidean_distance(q,X(i,:));
end
X = [X dist];
% Sort distances in ascending order
sort_X = sortrows(X,col+1,'ascend');
% Select k data points with the smallest distances 
ns = sort_X(1:k,1:end-1);
end