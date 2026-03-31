% =========================================================================
% 二元形式立方和分解的 Hankel 判别法 - 算法复现 (例 5.1)
% =========================================================================
clear; clc;

% 【例 5.1 数据录入】
% F(x,y) = 3x^6 + 15x^5y + 54x^4y^2 + 119x^3y^3 + 198x^2y^4 + 195xy^5 + 129y^6
d = 3; % 基底多项式次数，整体次数为 3d = 6
%A = [1,9,30,90,204,396,650,774,771,513]; % A_0 到 A_d+1
A=[1,2,3,4,5,6,7,8,9,10]
fprintf('=== 广义三阶(立方和) Waring 分解算法启动 ===\n');
fprintf('目标多项式次数 3d = %d\n', 3*d);
fprintf('输入原始系数 A_n: %s\n\n', mat2str(A));

% --- 第 1, 2 步：计算组合数与系数归一化 ---
a = normalize_coefficients_cubic(A, d);
fprintf('[第 1, 2 步] 提取并归一化系数 a_n:\n');
disp(round(a, 4));

% --- 第 3 步：构造 Hankel 矩阵并求秩 ---
try
    [r, c] = find_waring_rank_cubic(a, d);
    fprintf('[第 3 步] 成功判定广义 Waring 秩 r = %d\n', r);
    fprintf('          求解出递推方程系数 C = %s\n\n', mat2str(c', 4));
catch ME
    fprintf('算法终止：%s\n', ME.message);
    return;
end

% --- 第 4 步：求解特征方程，寻找基底参数 beta ---
beta = find_beta_params_cubic(c);
fprintf('[第 4 步] 特征多项式求根完毕，基底参数 beta_k 为:\n');
disp(beta.');

% --- 第 5 步：求解范德蒙系统，获取权重 lambda ---
lambda = solve_lambda_weights_cubic(a, beta, r);
fprintf('\n[第 5 步] 范德蒙系统求解完毕，对应权重 lambda_k 为:\n');
disp(lambda.');

% --- 第 6 步：结果展示 ---
fprintf('\n=== 最终立方和分解结果 ===\n');
for k = 1:r
    fprintf('项 %d: 权重 lambda = %.4f, 基底参数 beta = %.4f\n', k, lambda(k), beta(k));
    fprintf('      对应的基底 Q_%d(x,y) = x^%d + (%.4f)xy + (%.4f)y^%d\n', ...
            k, d, beta(k), beta(k)^2, d);
end





% --- 第 7 步：终极展开验算 ---
fprintf('\n=== 终极展开验算 ( \x03A3 \x03BB_i * Q_i^3 - F(x,y) ) ===\n');

% 1. 构造并输出基底 Q_k 的系数向量
% Q_k 对应的系数为 [1, beta_k, beta_k^2, ..., beta_k^d]
Q_coeffs = zeros(r, d + 1);
for k = 1:r
    for j = 0:d
        Q_coeffs(k, j+1) = beta(k)^j;
    end
    fprintf('Q_%d 提取的系数向量: %s\n', k, mat2str(round(Q_coeffs(k, :), 4)));
end

% 2. 利用卷积计算 sum(lambda_k * Q_k^3)
F_reconstructed = zeros(1, 3*d + 1); % 初始化重构的 F(x,y) 系数向量
for k = 1:r
    Q_k = Q_coeffs(k, :);
    
    % 多项式乘法在 MATLAB 中等价于系数向量的卷积
    Q_k_sq = conv(Q_k, Q_k);         % Q_k^2
    Q_k_cube = conv(Q_k_sq, Q_k);    % Q_k^3
    
    % 累加到总的多项式中
    F_reconstructed = F_reconstructed + lambda(k) * Q_k_cube;
end

% 3. 对比结果并计算残差 (Residual)
fprintf('\n重构的 \x03A3 \x03BB_i * Q_i^3 系数: %s\n', mat2str(round(F_reconstructed, 4)));
fprintf('原始的输入 F(x,y) 系数: %s\n', mat2str(round(A, 4)));

% 矩阵相减得到误差
residual = F_reconstructed - A;
fprintf('\n残差向量 (重构系数 - 原始系数): %s\n', mat2str(residual, 4));
fprintf('最大绝对误差 (Max Absolute Error): %e\n', max(abs(residual)));

if max(abs(residual)) < 1e-6
    fprintf('\n>>> 验证成功！重构误差极小，完美复原了原多项式！ <<<\n');
else
    fprintf('\n>>> 警告：存在较大的重构误差，请检查数值稳定性。 <<<\n');
end





