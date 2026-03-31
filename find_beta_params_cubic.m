function beta = find_beta_params_cubic(c)
    % 构造首一特征多项式 T(lambda) = lambda^r + c_{r-1}lambda^{r-1} + ... + c_0
    poly_coeffs = [1; flipud(c)];
    beta = roots(poly_coeffs);
    
    if length(unique(round(beta, 5))) < length(beta)
        warning('特征多项式存在重根，分解可能退化。');
    end
end