function intensity = model_3d_skewed(p, xyz)
    % p_distorted = [A, mu_x..z, sig_x..z, rho_xy..yz, alpha_x..z]
    p_distorted = p(1:10);
    alpha = p(11:13);
    mu = p(2:4);
    
    % Get the value of the un-skewed distorted Gaussian
    g_val = snap_helpers.model_3d_distorted(p_distorted, xyz);
    
    % Calculate the skew term (CDF)
    x_centered = xyz - mu;
    skew_term = 0.5 * (1 + erf((x_centered * alpha') / sqrt(2)));
    
    intensity = 2 * g_val .* skew_term;
end
