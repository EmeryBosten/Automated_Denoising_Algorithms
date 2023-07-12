function [dist] = euclidean_distance(point1,point2)
% calculating the Euclidean distance between two points in an n-dimensional space
% point1 and point2 are tuples with n coordinates
diff = point2 - point1;
sum_of_squares = (sum(diff.^2));
dist = sqrt(sum_of_squares);
end