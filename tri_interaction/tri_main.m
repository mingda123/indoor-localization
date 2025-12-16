% 基于三角关系的定位技术 - 多算法性能比较与CRLB分析
clear; clc; close all;
C = 299792458;
%% 参数设置
% 1. 设定三站坐标
sta_pos = [
    0, 0;
    1000, 0;
    500, 866
];

% 2. 设定移动台真实坐标
true_pos = [300; 400]; % 2x1向量

% 3. 计算真实距离
num_stations = size(sta_pos, 1);
true_distances = zeros(num_stations, 1);
for i = 1:num_stations
    true_distances(i) = norm(sta_pos(i,:)' - true_pos);
end

% 仿真参数
num_monte_carlo = 5000;  % Monte Carlo测试次数

% 噪声水平范围 (cm)
noise_levels_cm = logspace(0, 2, 10); % 10^0到 10^2 cm
% noise_levels_cm = linspace(1,100,100);
noise_levels_m = noise_levels_cm / 100; % 转换为米

% 要测试的算法列表
algorithms = {
    'LOP', ...
    'Gauss-Newton',...
    'ChanTOA'
};

% 初始化结果存储
rmse_results = zeros(length(algorithms), length(noise_levels_m));
mean_error_results = zeros(length(algorithms), length(noise_levels_m));
crlb_values = zeros(1, length(noise_levels_m));

%% 计算CRLB
fprintf('计算克拉美罗下界(CRLB)...\n');
for noise_idx = 1:length(noise_levels_m)
    noise_std_m = noise_levels_m(noise_idx);
    crlb_values(noise_idx) = calculate_CRLB(sta_pos, true_pos, noise_std_m);
    fprintf('噪声水平 %.1f cm: CRLB = %.4f m\n', noise_levels_cm(noise_idx), crlb_values(noise_idx));
end

%% 主测试循环 - 不同噪声水平
fprintf('\n开始多算法性能比较测试\n');

for noise_idx = 1:length(noise_levels_m)
    noise_std_m = noise_levels_m(noise_idx);
    fprintf('测试噪声水平: %.1f cm\n', noise_levels_cm(noise_idx));
    
    for algo_idx = 1:length(algorithms)
        algorithm = algorithms{algo_idx};
        errors = zeros(num_monte_carlo, 1);
        
        % Monte Carlo 测试
        for mc_idx = 1:num_monte_carlo
            % 产生带噪声的观测距离
            noise = noise_std_m * randn(num_stations, 1);
            measured_distances = true_distances + noise;
            
            % 使用不同算法估计位置
            switch algorithm
                case 'LOP'
                    estimated_pos = LOP(measured_distances, sta_pos);
                case 'Gauss-Newton'
                    estimated_pos = GaussNewtonMethod(measured_distances, sta_pos, 100, 1e-6);
                case 'ChanTOA'
                    estimated_pos = ChanTOA(sta_pos,measured_distances/C);
                otherwise
                    error('未知算法: %s', algorithm);
            end
            
            % 计算误差
            errors(mc_idx) = norm(estimated_pos - true_pos);
        end
        
        % 计算RMSE和平均误差
        rmse_results(algo_idx, noise_idx) = sqrt(mean(errors.^2));
        mean_error_results(algo_idx, noise_idx) = mean(errors);
        
        fprintf('  %s: RMSE = %.4f m (CRLB = %.4f m)\n', ...
            algorithm, rmse_results(algo_idx, noise_idx), crlb_values(noise_idx));
    end
    fprintf('\n');
end

%% 绘制性能比较图
figure('Position', [100, 100, 1000, 700]);

% 转换为厘米用于绘图
rmse_results_cm = rmse_results * 100;
crlb_values_cm = crlb_values * 100;

% 定义颜色和线型
colors = lines(length(algorithms) + 1);
line_styles = {'-', '--', ':', '-.'};
markers = {'o', 's', '^', 'd'};

% 绘制RMSE曲线
for algo_idx = 1:length(algorithms)
    loglog(noise_levels_cm, rmse_results_cm(algo_idx, :), ...
           'Color', colors(algo_idx, :), ...
           'LineStyle', line_styles{algo_idx}, ...
           'Marker', markers{algo_idx}, ...
           'MarkerSize', 6, ...
           'LineWidth', 2, ...
           'DisplayName', algorithms{algo_idx});
    hold on;
end

% 添加CRLB参考线 (理论下界)
loglog(noise_levels_cm, crlb_values_cm, 'k--', ...
       'LineWidth', 3, ...
       'DisplayName', 'CRLB理论下界');

% 图表美化
xlabel('测量噪声标准差 / cm', 'FontSize', 12, 'FontName', 'SimHei');
ylabel('RMSE / cm', 'FontSize', 12, 'FontName', 'SimHei');
title('TOA定位算法性能比较与CRLB分析', 'FontSize', 14, 'FontName', 'SimHei');
legend('Location', 'northwest', 'FontSize', 10, 'FontName', 'SimHei');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([10, 100]);
ylim([10, 1000]);

%% 性能分析表格
fprintf('\n=== 算法性能总结 (单位: cm) ===\n');
fprintf('噪声水平 ');
for algo_idx = 1:length(algorithms)
    fprintf('%14s ', algorithms{algo_idx});
end
fprintf('%14s\n', 'CRLB下界');

for noise_idx = 1:length(noise_levels_cm)
    fprintf('%8.1f ', noise_levels_cm(noise_idx));
    for algo_idx = 1:length(algorithms)
        fprintf('%14.2f ', rmse_results_cm(algo_idx, noise_idx));
    end
    fprintf('%14.2f\n', crlb_values_cm(noise_idx));
end

%% 效率分析
fprintf('\n=== 效率分析 ===\n');
for algo_idx = 1:length(algorithms)
    efficiency = crlb_values_cm ./ rmse_results_cm(algo_idx, :);
    avg_efficiency = mean(efficiency);
    fprintf('%s算法平均效率: %.2f%% (相对于CRLB)\n', ...
        algorithms{algo_idx}, avg_efficiency * 100);
end