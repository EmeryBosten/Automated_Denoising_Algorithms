function [folds] = split(data,k)
%function that splits data into k folds (of more or less equal lengths)
%data = matrix with x and y for n observations (ideally already shuffled)
%k = integer
%folds is a matrix which has 3 columns: x,y and an integer ranging from 1
%to k assigning fold
[row,~]=size(data);
vec = zeros([row,1]);
if rem(row,k) == 0
    lengthFolds = row/k;
    start = 0;
    for i=1:k
        vec((start*lengthFolds)+1:(start+1)*lengthFolds) = i;
        start = start + 1;
    end
else 
    r = rem(row,k);
    lengthFolds = floor(row/k);
    start = 0;
    for i=1:k
        vec(start*lengthFolds+1:(start+1)*lengthFolds) = i;
        start = start + 1;
    end
    vec(start*lengthFolds+1:end)=linspace(1,r,r);
end
folds = [data vec];
end