function [target_pos] = GaussNewtonMethod(D,sta_pos,max_iteration,tolerance)
% target_pos:目标估计位置 2*1
% D: 测站距离观测值n*1
% sta_pos: 测站位置 n*2
 

valid_idx = find(~isnan(D) & D>0);
n = length(valid_idx);

if(n < 3)
    warning("基站少于3台，无法进行定位");
    target_pos = [NaN;NaN];
    return;
end

% 设定初始值
target_pos = mean(sta_pos(valid_idx, :), 1)';
delta_pos = [inf;inf];
iteration = 0;

while norm(delta_pos) > tolerance && iteration < max_iteration
    H = zeros(n,2);
    Z = zeros(n,1);
    for k =1:n
        i = valid_idx(k);
        d = sqrt((target_pos(1) - sta_pos(i,1))^2+(target_pos(2) - sta_pos(i,2))^2);

        % 避免除以零
        if d > 0
            H(i, 1) = (target_pos(1) - sta_pos(i, 1)) / d;
            H(i, 2) = (target_pos(2) - sta_pos(i, 2)) / d;
        else
            H(i, 1) = 0;
            H(i, 2) = 0;
        end

        Z(i) = D(i) - d; 
    end

    delta_pos =(H'*H)\H'*Z;
    target_pos = target_pos + delta_pos;
    iteration = iteration + 1;
end
end