function [shuffledData] = shuffle(data)
%function that randomly shuffles input data
%data = matrix storing x and y for n observations 
[row,~] = size(data);
r = randperm(row);
shuffledData = zeros(size(data));
for k = 1:row
    shuffledData(k,:) = data(r(k),:);
end
end