% 三站定位算法 - 基于解析解法
function estimated_pos = LOP(measured_distances, sta_pos)
    
    % 提取三个基站的坐标
    x1 = sta_pos(1,1); y1 = sta_pos(1,2);
    x2 = sta_pos(2,1); y2 = sta_pos(2,2);
    x3 = sta_pos(3,1); y3 = sta_pos(3,2);
    
    % 提取三个距离测量值
    d1 = measured_distances(1);
    d2 = measured_distances(2);
    d3 = measured_distances(3);
    
    % 计算 k_i^2 = x_i^2 + y_i^2
    k2_sq = x2^2 + y2^2;
    k3_sq = x3^2 + y3^2;
    
    % 构建矩阵 H 和向量 b
    H = [x2, y2; 
         x3, y3];
    
    b = 0.5 * [k2_sq - d2^2 + d1^2;
               k3_sq - d3^2 + d1^2];
    
    % 求解线性方程组 H * x = b
    % 使用最小二乘解，以防矩阵接近奇异
    estimated_pos = H \ b;
    
    % 另一种实现方式：显式求逆（与图中公式一致）
    % 检查矩阵是否可逆
    if abs(det(H)) > 1e-10
        estimated_pos = inv(H) * b;
    else
        % 如果矩阵接近奇异，使用伪逆
        estimated_pos = pinv(H) * b;
    end
end