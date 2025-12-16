clc;clear;
%% 1.原始数据读取(两个静态点，一个动态)
static_file1 = "data/20251030_155850RTLS_log.txt";
static_file2 = "data/20251030_160011RTLS_log.txt";
dy_file = "data/20251030_160145RTLS_log.txt";

[static_anc1,static_epoch1] = UWB_fileread(static_file1);
[static_anc2,static_epoch2] = UWB_fileread(static_file2);
[dy_anc,dy_epoch] = UWB_fileread(dy_file);

%% 2.静态数据分析(二维x,y)
TS = [0,4.80;2.95,5.23];    % 静态参考真值

% 提取改正距离观测值并使用自编算法进行定位
MY_static1 = UWB_postioning(static_anc1,static_epoch1);
fprintf('--- 开始处理静态点1 ---\n');
static_analysis(MY_static1, TS(1,:), static_epoch1, '静态点1 - XY平面分析', true);
MY_static2 = UWB_postioning(static_anc2,static_epoch2);
fprintf('--- 开始处理静态点2 ---\n');
static_analysis(MY_static2, TS(2,:), static_epoch2, '静态点2 - XY平面分析', true);

%% 3.动态数据分析
fprintf('\n--- 开始处理动态数据 ---\n');
MY_dynamic = UWB_postioning(dy_anc, dy_epoch);
% 调用新的动态数据分析函数
dynamic_analysis(MY_dynamic, dy_epoch, '动态轨迹量化分析', true);


%% 静态点误差分析
function static_analysis(position_results_3d, truth_pos_2d, epochs, figure_title, plot_le)
% 静态点误差分析

    if nargin < 5, plot_le = false; end
    
    % 数据投影
    my_results_2d = position_results_3d(:, 1:2);
    le_results_2d = [];
    if plot_le
        num_epochs = length(epochs);
        le_results_2d = NaN(num_epochs, 2);
        for i = 1:num_epochs
            if ~isempty(epochs(i).LE) && isfield(epochs(i).LE, 'x')
                le_results_2d(i, 1) = epochs(i).LE.x;
                le_results_2d(i, 2) = epochs(i).LE.y;
            end
        end
    end
    
    figure('Name', figure_title, 'NumberTitle', 'off', 'WindowState', 'maximized');
    % XY平面轨迹
    subplot(2, 2, 1);
    hold on; grid on;
    plot(my_results_2d(:,1), my_results_2d(:,2), 'b.-', 'DisplayName', '自解算轨迹');
    plot(truth_pos_2d(1), truth_pos_2d(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', '参考真值');
    if plot_le, plot(le_results_2d(:,1), le_results_2d(:,2), 'm.-', 'DisplayName', 'UWB内置算法轨迹(未平滑)'); end
    title('定位成果图 (XY平面)'); xlabel('X (m)'); ylabel('Y (m)');
    legend('Location', 'best'); axis equal;

    % 平面总误差
    subplot(2, 2, 2);
    hold on; grid on;
    errors_my_algo_2d = sqrt(sum((my_results_2d - truth_pos_2d).^2, 2));
    plot(errors_my_algo_2d, 'b.-', 'DisplayName', '自解算总误差');
    if plot_le
        errors_le_2d = sqrt(sum((le_results_2d - truth_pos_2d).^2, 2));
        plot(errors_le_2d, 'm.-', 'DisplayName', 'UWB内置算法总误差(未平滑)');
    end
    title('平面定位总误差'); xlabel('历元序号'); ylabel('误差 (m)');
    legend('Location', 'best');

    % X方向和Y方向误差分解
    % X方向
    subplot(2, 2, 3);
    hold on; grid on;
    errors_my_algo_x = my_results_2d(:,1) - truth_pos_2d(1);
    plot(errors_my_algo_x, 'Color', [0 0.4470 0.7410], 'DisplayName', '自解算X误差'); % 蓝色
    plot([1, length(errors_my_algo_x)], [0, 0], 'k--', 'HandleVisibility', 'off'); % 零线
    if plot_le
        errors_le_x = le_results_2d(:,1) - truth_pos_2d(1);
        plot(errors_le_x, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'UWB内置算法X误差(未平滑)'); % 橙色
    end
    title('X方向误差 (测量值 - 真值)'); xlabel('历元序号'); ylabel('误差 (m)');
    legend('Location', 'best');

    % Y方向
    subplot(2, 2, 4);
    hold on; grid on;
    errors_my_algo_y = my_results_2d(:,2) - truth_pos_2d(2);
    plot(errors_my_algo_y, 'Color', [0 0.4470 0.7410], 'DisplayName', '自解算Y误差');
    plot([1, length(errors_my_algo_y)], [0, 0], 'k--', 'HandleVisibility', 'off');
    if plot_le
        errors_le_y = le_results_2d(:,2) - truth_pos_2d(2);
        plot(errors_le_y, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'UWB内置算法Y误差(未平滑)');
    end
    title('Y方向误差 (测量值 - 真值)'); xlabel('历元序号'); ylabel('误差 (m)');
    legend('Location', 'best');

    fprintf('--- XY平面详细分析报告: %s ---\n', figure_title);
    fprintf('自解算算法 (XY平面):\n');
    fprintf('  总误差(RMSE): %.3f m\n', sqrt(mean(errors_my_algo_2d.^2, 'omitnan'))); % RMSE
    fprintf('  X方向误差: 均值(Bias)=%.3f m, 标准差(Precision)=%.3f m\n', mean(errors_my_algo_x, 'omitnan'), std(errors_my_algo_x, 'omitnan'));
    fprintf('  Y方向误差: 均值(Bias)=%.3f m, 标准差(Precision)=%.3f m\n', mean(errors_my_algo_y, 'omitnan'), std(errors_my_algo_y, 'omitnan'));
    if plot_le
        fprintf('UBW内置算法(未平滑)(XY平面):\n');
        fprintf('  总误差(RMSE): %.3f m\n', sqrt(mean(errors_le_2d.^2, 'omitnan')));
        fprintf('  X方向误差: 均值(Bias)=%.3f m, 标准差(Precision)=%.3f m\n', mean(errors_le_x, 'omitnan'), std(errors_le_x, 'omitnan'));
        fprintf('  Y方向误差: 均值(Bias)=%.3f m, 标准差(Precision)=%.3f m\n', mean(errors_le_y, 'omitnan'), std(errors_le_y, 'omitnan'));
    end
    fprintf('--- 分析结束 ---\n');
end
%% 动态误差分析
function dynamic_analysis(position_results_3d, epochs, figure_title, plot_le)
% 动态误差分析

    if nargin < 4, plot_le = false; end

    % 数据投影
    my_results_2d = position_results_3d(:, 1:2);
    le_results_2d = [];
    if plot_le
        num_epochs = length(epochs);
        le_results_2d = NaN(num_epochs, 2);
        for i = 1:num_epochs
            if ~isempty(epochs(i).LE) && isfield(epochs(i).LE, 'x')
                le_results_2d(i, 1) = epochs(i).LE.x;
                le_results_2d(i, 2) = epochs(i).LE.y;
            end
        end
    end

    figure('Name', figure_title, 'NumberTitle', 'off', 'WindowState', 'maximized');
    % 动态轨迹
    hold on; grid on;
    plot(my_results_2d(:,1), my_results_2d(:,2), 'b.-', 'DisplayName', '自解算动态轨迹');
    if plot_le, plot(le_results_2d(:,1), le_results_2d(:,2), 'm.-', 'DisplayName', 'UWB内置算法动态轨迹(未平滑)'); end
    title('动态定位轨迹对比 (XY平面)'); xlabel('X (m)'); ylabel('Y (m)');
    legend; axis equal;
  
    fprintf('--- 动态轨迹分析报告: %s ---\n', figure_title);
    fprintf('--- 分析结束 ---\n');
end