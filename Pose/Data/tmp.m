% 读取并显示图像
img = imread('pic.jpg');
imshow(img);

% 获取用户点击的坐标点
n = 4;
[x, y] = ginput(n);  % n为要获取的点数

% 显示坐标
disp('获取的坐标点:');
for i = 1:length(x)
    fprintf('点%d: (%.1f, %.1f)\n', i, x(i), y(i));
end

% 在图上标记选中的点
hold on;
plot(x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
hold off;