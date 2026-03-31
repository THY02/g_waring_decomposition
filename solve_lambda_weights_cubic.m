function lambda = solve_lambda_weights_cubic(a, beta, r)
    V = zeros(r, r);
    for i = 1:r
        for j = 1:r
            V(i, j) = beta(j)^(i-1); 
        end
    end
    a_sub = a(1:r)';
    lambda = V \ a_sub;
end