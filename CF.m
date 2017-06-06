function [ ret ] = CF( x, f, d )
%CF Summary of this function goes here
%   Detailed explanation goes here
    % x : 3 x 1
    % f : 3 x number of feature
    ret = d(:,f(1,:)==x(1,1) & f(2,:)==x(2,1));
    disp(ret);
end

