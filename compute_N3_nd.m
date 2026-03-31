function N = compute_N3_nd(n, d)
    % COMPUTE_N3_ND 计算等比基底立方展开的组合系数 N_3(n,d)
    % 等价于求 (1 + t + ... + t^d)^3 展开式中 t^n 的系数
    
    if n < 0 || n > 3*d
        N = 0;
        return;
    end
    
    % 构造多项式 (1 + t + ... + t^d) 的系数向量 [1, 1, ..., 1]
    poly_base = ones(1, d + 1);
    
    % 利用 MATLAB 的 conv 函数进行卷积，求三次方展开的系数
    poly_sq = conv(poly_base, poly_base);      % 平方
    poly_cube = conv(poly_sq, poly_base);      % 立方
    
    % 返回对应 t^n 的系数（MATLAB 索引从 1 开始，所以取 n+1）
    N = poly_cube(n + 1);
end