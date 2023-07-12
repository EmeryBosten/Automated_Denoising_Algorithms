function [V] = modoneNNestimator(Y,X,type)
% computes modified 1-nearest neighbor residual variance estimator
% Implementation follows "Residual variance estimation in machine learning"
% by Elia Liitiainen, Michel Verleysen, Francesco Corona, and Amaury Lendasse
% Author modified to median for more robust estimation in the presence of
% sudden changes in the signal by sharp peaks
k = 3;
len = length(Y);
Yn1 = zeros([len, 1]);
Yn2 = zeros([len, 1]);

for i = 1:len
    n = knn(X,i,k);
    n(1)=[];
    Yn1(i) = Y(n(1));
    Yn2(i) = Y(n(2));

end
diff1 = abs(Y - Yn1);
diff2 = abs(Y - Yn2);
mult_diff = diff1.*diff2;
if strcmp(type,'median')
    V = median(mult_diff); % take median because of peaks in chromatograms
% peaks are sudden changes in value and will act as outliers in mean 
% need more robust estimation 
elseif strcmp(type,'mean')
    sum_mult_diff = sum(mult_diff);
    V = sum_mult_diff/(len);
end
end