function a = normalize_coefficients_cubic(A, d)
    % NORMALIZE_COEFFICIENTS_CUBIC 将原始系数 A_n 归一化为 a_n
    % 输入: A 为 [A_0, A_1, ..., A_{3d}]
    
    a = zeros(size(A));
    for n = 0:(3*d)
        N_val = compute_N3_nd(n, d);
        a(n+1) = A(n+1) / N_val; 
    end
end