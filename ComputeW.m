function [ W ] = ComputeW( Intensity, X, X_prime, X_tilde_prime, Hn_1, f_x, d_x, f_xp, d_xp)
%COMPUTEW Summary of this function goes here
%   Detailed explanation goes here
    % X = (3 x n) matrix
    % X_prime = (3 x n) matrix
    [r, c] = size(X);
    W = zeros(c,c);
    [r_i, c_i] = size(Intensity);
    L_I = reshape(Intensity, 1, r_i * c_i);
    p_10_i = 1 / (max(L_I) - min(L_I));
    for i=1:c
        % P(xi, z = 0 | Hn_1) = 1/U
        % P(xi | z = 1, Hn_1) = PDF( x, x_prime, H, sigma )
        x_i = X(:,i);
        x_p = X_prime(:,i);
        x_tilde_p = X_tilde_prime(:,i);
        
        temp2 = H_n1\x_tilde_p;
        temp4 = (x_i - temp2).' * (x_i - temp2);
        
        p_4 = PDF(x_i, x_p, H_n1, 1.0, f_x, d_x, f_xp, d_xp);
        p_17 = exp(temp4);
        W(i,i) = p_17 / (p_10_i + p_4);
    end
end

