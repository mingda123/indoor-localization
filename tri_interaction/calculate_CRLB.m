function crlb = calculate_CRLB(sta_pos, true_pos, noise_std)
 
    
    num_stations = size(sta_pos, 1);
    true_pos = true_pos(:); % 确保是列向量
    
    % 计算真实距离和方向向量
    true_distances = zeros(num_stations, 1);
    direction_vectors = zeros(num_stations, 2);
    
    for i = 1:num_stations
        diff = sta_pos(i,:)' - true_pos;
        true_distances(i) = norm(diff);
        if true_distances(i) > 0
            direction_vectors(i,:) = (diff / true_distances(i))';
        end
    end
    
    % 构建Fisher信息矩阵
    FIM = zeros(2,2);
    for i = 1:num_stations
        g = direction_vectors(i,:)';
        FIM = FIM + (g * g') / (noise_std^2);
    end
    
    % CRLB是Fisher信息矩阵的逆
    if cond(FIM) > 1e10
        crlb_matrix = pinv(FIM); % 使用伪逆避免数值问题
    else
        crlb_matrix = inv(FIM);
    end
    
    % CRLB是位置估计误差方差的下界
    crlb = sqrt(trace(crlb_matrix)); % RMSE形式
end