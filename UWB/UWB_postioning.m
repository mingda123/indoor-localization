function [target_pos] = UWB_postioning(anc,epoch) 
% 解析UWB数据进行定位 
% 输入: 
%   anc - 基站信息 
%   epoch - 每历元观测定位数据
% 输出: 
%   target_pos - 每历元自定位结果 

    % 获取基站坐标（x,y） 
    n = size(anc,2); 
    sta_pos = zeros(n,3); 
    for i = 1:n 
        id = anc(i).id; 
        sta_pos(id+1,1) = anc(i).x; 
        sta_pos(id+1,2) = anc(i).y; 
        sta_pos(id+1,3) = anc(i).z; 
    end 
    len = size(epoch,2); 
    target_pos = zeros(len,3);
    for i =1:len 
        D = NaN(n,1); 
        for j =1:size(epoch(i).RR,2) 
            range_meas = epoch(i).RR(j).corrected_range / 1000;
            if(~isempty(range_meas))
                D(j,1) = range_meas; 
            end 
        end 
        target_pos(i,:) = GaussNewton_Refined(D,sta_pos,100,1e-6);
    end 
end

function [initial_pos] = LSM_initial_pos(D, sta_pos)
% 初值求解
% D: 测站距离观测值 n*1
% sta_pos: 测站位置 n*3
% initial_pos: 目标估计位置的初始解 3*1
    n = size(sta_pos, 1); % 获取基站数量

    % 至少需要3个基站才能进行2D定位/4基站3D定位
    if n < 4
        error('至少需要4个基站才能计算初始位置。');
    end

    % 初始化A、b
    A = zeros(n, 4);
    b = zeros(n, 1);
    for i = 1:n
        Xi = sta_pos(i, 1);
        Yi = sta_pos(i, 2);
        Zi = sta_pos(i, 3);
        Di = D(i);
        A(i, :) = [-2*Xi, -2*Yi, -2*Zi, 1];
        b(i) = Di^2 - Xi^2 - Yi^2 - Zi^2;
    end
    x_solved = (A'*A)^-1 * A'*b;
    initial_pos = [x_solved(1); x_solved(2);x_solved(3)]; % 只取前两个元素
end

function [target_pos] = GaussNewton_Refined(D, sta_pos, max_iteration, tolerance)
% 迭代计算
% D: 测站距离观测值 n*1
% sta_pos: 测站位置 n*3
% max_iteration: 最大迭代次数
% tolerance: 收敛阈值
% target_pos: 最终的目标估计位置 3*1

    valid_idx = find(~isnan(D) & D>0);
    n = length(valid_idx);
    
    if(n < 4)
        warning("基站少于4台，无法进行定位");
        target_pos = [NaN;NaN;NaN];
        return;
    end
    
    % 设定初始值
    target_pos = LSM_initial_pos(D,sta_pos);
    delta_pos = [inf;inf;inf];
    iteration = 0;
    
    while norm(delta_pos) > tolerance && iteration < max_iteration
        H = zeros(n,3);
        Z = zeros(n,1);
        for k =1:n
            i = valid_idx(k);
            d = sqrt((target_pos(1) - sta_pos(i,1))^2+(target_pos(2) - sta_pos(i,2))^2 + (target_pos(3) - sta_pos(i,3))^2);
    
            % 避免除以零
            if d > 0
                H(i, 1) = (target_pos(1) - sta_pos(i, 1)) / d;
                H(i, 2) = (target_pos(2) - sta_pos(i, 2)) / d;
                H(i, 3) = (target_pos(3) - sta_pos(i, 3)) / d;
            else
                H(i, :) = 0;
            end
    
            Z(i) = D(i) - d; 
        end
    
        delta_pos =(H'*H)\H'*Z;
        target_pos = target_pos + delta_pos;
        iteration = iteration + 1;
    end
end