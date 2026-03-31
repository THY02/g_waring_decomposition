function [r, c] = find_waring_rank_cubic(a, d)
    % 适用范围：2r - 1 <= 3d
    max_r = floor((3*d + 1) / 2); 
    
    for r = 1:max_r
        % 1. 构造 r x r 的 Hankel 矩阵 H_r
        H_r = zeros(r, r);
        for i = 1:r
            for j = 1:r
                H_r(i, j) = a((i-1) + (j-1) + 1);
            end
        end
        
        % 2. 构造等式右边的向量 -A_r
        A_r = zeros(r, 1);
        for i = 1:r
            A_r(i) = a((i-1) + r + 1);
        end
        
        % 3. 检查矩阵是否接近奇异（不可逆）
        if cond(H_r) > 1e10
            continue; 
        end
        
        % 4. 求解线性方程组 H_r * c = -A_r
        c = H_r \ (-A_r);
        
        % 5. 剩余递推式全局检验 a_{n+r} + c_{r-1}a_{n+r-1} + ... + c_0 a_n = 0
        valid = true;
        for n = 0:(3*d - r)
            sum_val = a(n + r + 1);
            for k = 0:(r-1)
                sum_val = sum_val + c(k+1) * a(n + k + 1);
            end
            if abs(sum_val) > 1e-6
                valid = false;
                break;
            end
        end
        
        if valid
            return;
        end
    end
    error('在低秩假设 (2r-1 <= 3d) 下未找到解，请考虑使用催化矩阵方法。');
end