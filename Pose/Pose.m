clear; clc; close all;

%% 1. 输入数据

% % 相机内参
% f = 3018.65;
% x0 = 1508.85;
% y0 = 1953.68;
% K = [f, 0, x0;
%      0, f, y0;
%      0, 0, 1];
% k1 = 0.0102086;k2 = -0.0810434;k3 = 0.0886927;
% p1 = 0.000484137;p2 = -3.54084e-005;
% 
% % 4个特征点的像素坐标 (Upix, Vpix)
% imagePoints = [1077, 3609;
%                693,  736;
%                2189, 914;
%                2399, 2935];
% 
% % 对应的4个特征点的世界坐标 (Xw, Yw, Zw)
% objectPointsWorld = [-3.3, 121.5, 61.4;
%                      -3.3, 121.5, 184.7;
%                      90,   121.5, 184.7;
%                      90,   121.5, 61.4];
f = (2809.27724 + 2714.81215)/2;
x0 = 4096/2+(-83.2887-41.474)/2;
y0 = 3072/2+(18.9526+5.74873)/2;
K = [f, 0, x0;
     0, f, y0;
     0, 0, 1];
k1 = (-0.0161696-0.0224794)/2;
k2 = (0.17964+0.143419)/2;
k3 = (-0.217107-0.152107)/2;
p1 = (0.00238636+0.000701319)/2;
p2 = (0.000416543+0.00067288)/2;
imagePoints = [1395, 767;
               3163,  851;
               1347, 1999;
               3207, 2007];
objectPointsWorld = [14.35, 161.85, 112.48;
                     14.35+29.98, 161.85, 112.48;
                     14.35, 161.85, 92.45;
                     14.35+29.98, 161.85, 92.45];

fprintf('--- 输入数据 ---\n');
disp('相机内参矩阵 K:'); disp(K);
disp('像点坐标 (imagePoints):'); disp(imagePoints);
disp('物方世界坐标 (objectPointsWorld):'); disp(objectPointsWorld);



%% 2. 建立物方平面坐标系，使得四角坐标Z=0
meanPointsWorld = mean(objectPointsWorld,1);
A_centered = objectPointsWorld - meanPointsWorld;
[~,~,V_svd] = svd(A_centered'*A_centered);
Rpl = V_svd;
tpl = meanPointsWorld';
if det(Rpl) < 0
    fprintf('检测到Rpl是镜像变换，翻转第三列使其成为旋转矩阵\n');
    Rpl(:,3) = -Rpl(:,3);  % 翻转法向量方向
end
fprintf('修正后Rpl行列式: %.6f\n', det(Rpl));

Xw = objectPointsWorld'; 
Xpl_h = Rpl' * (Xw - tpl);
Xpl = Xpl_h(1:2, :); % 只取前两行(X, Y)作为平面坐标
%% 3. 像素坐标转为归一化坐标系下像点坐标
Xpix = [imagePoints'; ones(1, size(imagePoints,1))]; 

Xnp = K \ Xpix;
%% 4. 根据平面坐标系坐标与归一化坐标系下像点坐标，估计平面坐标系下Hpl矩阵
A = zeros(8,9);

for i = 1:4
    Xpli = Xpl(1, i); Ypli = Xpl(2, i);
    Xnpi = Xnp(1, i); Ynpi = Xnp(2, i);
    
    A(2*i-1, :) = [Xpli, Ypli, 1, 0, 0, 0, -Xpli*Xnpi, -Ypli*Xnpi, -Xnpi];
    A(2*i, :)   = [0, 0, 0, Xpli, Ypli, 1, -Xpli*Ynpi, -Ypli*Ynpi, -Ynpi];
end

% 使用SVD求解 Ah=0
[~, ~, V] = svd(A);
h = V(:, end); 
Hpl = reshape(h, 3, 3)';
%% 5. 分解Hpl，得到Rh和th
scale = (norm(Hpl(:,1)) + norm(Hpl(:,2))) / 2;
r1_raw = Hpl(:,1) / norm(Hpl(:,1));
r2_raw = Hpl(:,2) / norm(Hpl(:,2));
th = Hpl(:,3) / scale;

r1 = r1_raw;
r2 = r2_raw - dot(r2_raw, r1) * r1;
r2 = r2 / norm(r2);
r3 = cross(r1, r2);

Rh = [r1, r2, r3];
if det(Rh) < 0
    fprintf('检测到镜像变换，翻转第三列\n');
    Rh(:,3) = -Rh(:,3);  % 翻转r3方向
end

fprintf('修正后Rh行列式: %.6f\n', det(Rh));
fprintf('修正后Rh正交性: %.6e\n', norm(Rh*Rh' - eye(3), 'fro'));
%% 6. 将Rh和th转换到世界坐标系下
% 旋转矩阵 R (世界 -> 相机)
R = Rh * Rpl';

% 位移向量 t (世界 -> 相机)
t = th - R * tpl;

% 相机中心位置 C (在世界坐标系下)
C = -R' * t;
%% 7. 输出初值结果
fprintf('\n--- 初值计算结果 ---\n');
% 将旋转矩阵转换为欧拉角 (Omega, Phi, Kappa)
phi_rad = asin(-R(3,1));
omega_rad = atan2(R(3,2)/cos(phi_rad), R(3,3)/cos(phi_rad));
kappa_rad = atan2(R(2,1)/cos(phi_rad), R(1,1)/cos(phi_rad));
opk_deg_init = rad2deg([omega_rad, phi_rad, kappa_rad]);

fprintf('相机位置 XYZ (初始值/cm): [%.6f, %.6f, %.6f]\n', C(1), C(2), C(3));
fprintf('相机姿态 OPK (初始值/deg): [%.6f, %.6f, %.6f]\n', opk_deg_init);

%% 8. 迭代优化
% 设置初值
Xs = C(1); Ys = C(2); Zs = C(3);
% 保持初始的角元素值
omega = omega_rad; phi = phi_rad; kappa = kappa_rad;

% 迭代设置
max_iter = 50;
convergence_threshold = 1e-6;
eps_ang = 1e-7;               % 数值差分步长

fprintf('\n--- 开始迭代优化 ---\n');
fprintf('%-5s %-12s %-12s\n', '迭代', '改正数范数', '误差范数');
fprintf('-----------------------------------\n');

% 由 OPK 生成旋转矩阵（按 R = Rz(kappa) * Ry(phi) * Rx(omega)）
rot_from_opk = @(om, ph, ka) ...
    ( [cos(ka)*cos(ph), cos(ka)*sin(ph)*sin(om)-sin(ka)*cos(om), cos(ka)*sin(ph)*cos(om)+sin(ka)*sin(om);
       sin(ka)*cos(ph), sin(ka)*sin(ph)*sin(om)+cos(ka)*cos(om), sin(ka)*sin(ph)*cos(om)-cos(ka)*sin(om);
       -sin(ph),        cos(ph)*sin(om),                         cos(ph)*cos(om)] );


for iter = 1:max_iter
    % 计算旋转矩阵 R
    co = cos(omega); so = sin(omega);
    cp = cos(phi);   sp = sin(phi);
    ck = cos(kappa); sk = sin(kappa);
    
    R = rot_from_opk(omega, phi, kappa);

    % 初始化
    num_points = size(objectPointsWorld, 1);
    A = zeros(2 * num_points, 6);
    L = zeros(2 * num_points, 1);

    for i = 1:num_points
        % 当前处理的点的坐标
        X = objectPointsWorld(i, 1); 
        Y = objectPointsWorld(i, 2); 
        Z = objectPointsWorld(i, 3);
        x_obs = imagePoints(i, 1);    
        y_obs = imagePoints(i, 2);
        
        % 计算像空间辅助坐标
        v = [X - Xs; Y - Ys; Z - Zs];
        Xc = R * v;
        X_bar = Xc(1); Y_bar = Xc(2); Z_bar = Xc(3);
     
        % 计算近似像点坐标 x_calc, y_calc
        x_calc = x0 + f * X_bar / Z_bar;
        y_calc = y0 + f * Y_bar / Z_bar;
        
        % 误差向量计算
        L(2*i-1) = x_obs - x_calc;
        L(2*i)   = y_obs - y_calc;
        
        % 线元素偏导
        a11 = f * ( -R(1,1)*Z_bar - X_bar*(-R(3,1)) ) / (Z_bar^2);
        a12 = f * ( -R(1,2)*Z_bar - X_bar*(-R(3,2)) ) / (Z_bar^2);
        a13 = f * ( -R(1,3)*Z_bar - X_bar*(-R(3,3)) ) / (Z_bar^2);

        a21 = f * ( -R(2,1)*Z_bar - Y_bar*(-R(3,1)) ) / (Z_bar^2);
        a22 = f * ( -R(2,2)*Z_bar - Y_bar*(-R(3,2)) ) / (Z_bar^2);
        a23 = f * ( -R(2,3)*Z_bar - Y_bar*(-R(3,3)) ) / (Z_bar^2);
        
        % omega
        R_om = rot_from_opk(omega + eps_ang, phi, kappa);
        Xc_p = R_om * v;
        x_p = x0 + f * (Xc_p(1)/Xc_p(3));
        y_p = y0 + f * (Xc_p(2)/Xc_p(3));
        a14 = (x_p - x_calc) / eps_ang;
        a24 = (y_p - y_calc) / eps_ang;
        
        % phi
        R_ph = rot_from_opk(omega, phi + eps_ang, kappa);
        Xc_p = R_ph * v;
        x_p = x0 + f * (Xc_p(1)/Xc_p(3));
        y_p = y0 + f * (Xc_p(2)/Xc_p(3));
        a15 = (x_p - x_calc) / eps_ang;
        a25 = (y_p - y_calc) / eps_ang;
        
        % kappa
        R_ka = rot_from_opk(omega, phi, kappa + eps_ang);
        Xc_p = R_ka * v;
        x_p = x0 + f * (Xc_p(1)/Xc_p(3));
        y_p = y0 + f * (Xc_p(2)/Xc_p(3));
        a16 = (x_p - x_calc) / eps_ang;
        a26 = (y_p - y_calc) / eps_ang;
        
        % 填充A矩阵
        A(2*i-1, :) = [a11, a12, a13, a14, a15, a16];
        A(2*i, :)   = [a21, a22, a23, a24, a25, a26];
    end

    % 求取改正数dX
    damping = 1e-4;
    N = A' * A + damping * eye(6);
    dX = N \ (A' * L);
    
    % 更新参数
    Xs = Xs + dX(1);
    Ys = Ys + dX(2);
    Zs = Zs + dX(3);
    omega = omega + dX(4);
    phi   = phi   + dX(5);
    kappa = kappa + dX(6);

    % 显示迭代进程
    fprintf('%-5d %-12.4e %-12.4e\n', iter, norm(dX), norm(L));
    
    % 检查收敛性
    if norm(dX) < convergence_threshold
        fprintf('-----------------------------------\n');
        fprintf('成功收敛于第 %d 次迭代。\n', iter);
        break;
    end
end

if iter == max_iter && norm(dX) >= convergence_threshold
    fprintf('-----------------------------------\n');
    fprintf('警告：已达到最大迭代次数，但未收敛。\n');
end

% 最终优化后的参数
fprintf('\n--- 最终优化结果 ---\n');
fprintf('优化后相机位置 XYZ(cm): [%.6f, %.6f, %.6f]\n', Xs, Ys, Zs);
fprintf('优化后相机姿态 OPK (deg): [%.6f, %.6f, %.6f]\n', ...
        rad2deg(omega), rad2deg(phi), rad2deg(kappa));