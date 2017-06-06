function [ ret ] = PDF( x, x_prime, H, sigma, f_x, d_x, f_xp, d_xp )
%PDF Summary of this function goes here
%   Detailed explanation goes here
    temp = CF(x, f_x, d_x) - CF(H\x_prime, f_xp, d_xp);
    temp2 = (temp.' * temp) / (2*sigma*sigma);
    ret = exp(-temp2) / (sqrt(2*pi)*sigma);
end

