function corners = manual_select_corners(img)
    % 手动选择A4纸四个角点
    % 顺序：左上 → 右上 → 右下 → 左下
    
    figure;
    imshow(img);
    title('请按顺序点击A4纸的四个角点：左上 → 右上 → 右下 → 左下');
    
    corners = zeros(4, 2);
    
    for i = 1:4
        [x, y] = ginput(1);
        corners(i, 1) = x;
        corners(i, 2) = y;
        
        % 在图像上标记选中的点
        hold on;
        plot(x, y, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
        text(x, y, sprintf('  %d', i), 'Color', 'red', 'FontSize', 12);
        hold off;
    end
    
    % 连接角点显示四边形
    hold on;
    plot([corners(:,1); corners(1,1)], [corners(:,2); corners(1,2)], 'g-', 'LineWidth', 2);
    hold off;
end