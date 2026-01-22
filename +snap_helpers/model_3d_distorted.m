function intensity = model_3d_distorted(p, xyz)
    % p = [A, mu_x, mu_y, mu_z, sig_x, sig_y, sig_z, rho_xy, rho_xz, rho_yz]
    A = p(1);
    mu = p(2:4);
    sig = p(5:7);
    rho = p(8:10);
    
    % Construct covariance matrix
    Sigma = [sig(1)^2,              rho(1)*sig(1)*sig(2),   rho(2)*sig(1)*sig(3);
             rho(1)*sig(1)*sig(2),   sig(2)^2,               rho(3)*sig(2)*sig(3);
             rho(2)*sig(1)*sig(3),   rho(3)*sig(2)*sig(3),   sig(3)^2];
    
    % Ensure covariance matrix is positive semi-definite
    [~, flag] = chol(Sigma);
    if flag ~= 0
        intensity = zeros(size(xyz, 1), 1); % Return zero if matrix is invalid
        return;
    end
    
    x_centered = xyz - mu;
    
    intensity = A * exp(-0.5 * sum((x_centered / Sigma) .* x_centered, 2));
end
