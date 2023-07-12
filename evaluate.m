function [score] = evaluate(ztest,ytest,type)
%evaluate on test set
%ztest = interpolated estimation on test set
%ytest = hold back test set
% score = sum((ytest - ztest).^2)/length(ytest);
if strcmp(type,'median')
    score = median((ytest - ztest).^2);
elseif strcmp(type,'mean')
    score = sum((ytest - ztest).^2)/length(ytest);
end
end