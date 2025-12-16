function [target_pos] = ChanTOA(sta_pos,toa)
% Chan氏求解法
% target_pos:目标位置(2*1)
% sta_pos:测站位置(n*2)
% toa:信号从目标点到达各个基站的到达时间(n*1)

valid_idx = find(~isnan(toa) & toa>0);
n = length(valid_idx);
if(n < 3)
    warning("基站少于3台，无法进行定位");
    target_pos = [NaN;NaN];
    return;
end

% 初始化
C = 299792458;   % 光速

% 计算相应的矩阵
dr = zeros(n-1,1);
d = zeros(n,1);d(1) = sqrt(sta_pos(1,1)^2+sta_pos(1,2)^2);
D = zeros(n-1,1);
A = zeros(n-1,2);
Dr = zeros(n-1,1);
for i =2:n
    dr(i-1) = C*(toa(i) - toa(1));
    d(i) = sqrt(sta_pos(i,1)^2+sta_pos(i,2)^2);
    D(i-1) = (d(i)^2 - d(1)^2 - dr(i-1)^2)/2;

    A(i-1,1) = sta_pos(i,1) - sta_pos(1,1);
    A(i-1,2) = sta_pos(i,2) - sta_pos(1,2);

    Dr(i-1) = -dr(i-1);
end

% a1,a2
temp_a = (A'*A)^-1*A'* Dr;
a1 = temp_a(1);a2 = temp_a(2);

% b1,b2
temp_b = (A'*A)^-1*A'* D;
b1 = temp_b(1);b2 = temp_b(2);

% I,J,K,r
x1 = sta_pos(1,1);y1 = sta_pos(1,2);
I = a1^2+a2^2-1;
J = a1 * (b1 - x1) + a2 * (b2 - y1);
K = (x1 - b1)^2 + (y1 - b2)^2;
r = -(J+sqrt(J^2-I*K))/I;

% 计算目标点
target_pos = [a1 * r + b1;a2 * r + b2];
end