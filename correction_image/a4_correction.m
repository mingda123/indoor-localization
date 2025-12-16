function corrected_img = a4_correction_fixed(img, corners,dpi,width_cm,height_cm)
    % 修正版的A4纸图像校正
    
    % 定义校正空间
    width_pixels = width_cm * dpi;
    height_pixels = height_cm * dpi;
    fprintf('校正后图像尺寸: %d x %d\n', width_pixels, height_pixels);
    
    % 定义目标角点坐标（校正后的A4纸四个角）
    dst_corners = [1, 1; 
                   width_pixels, 1; 
                   width_pixels, height_pixels; 
                   1, height_pixels];
    
    % 计算H阵
    H = compute_homography(dst_corners, corners);
    
    % 检查H阵
    fprintf('单应性矩阵H:\n');
    disp(H);
    
    % 填充像素
    corrected_img = zeros(height_pixels, width_pixels, 3, 'uint8');
    
    % 获取原始图像尺寸
    [orig_height, orig_width, ~] = size(img);
    fprintf('原始图像尺寸: %d x %d\n', orig_width, orig_height);
    
    filled_pixels = 0;
    
    for y_corr = 1:height_pixels
        for x_corr = 1:width_pixels
            % 计算对应的原始图像坐标
            src_pt = H * [x_corr; y_corr; 1];
            src_pt = src_pt / src_pt(3);
            
            x_img = src_pt(1);
            y_img = src_pt(2);
            
            % 边界检查
            if x_img >= 1 && x_img <= orig_width && y_img >= 1 && y_img <= orig_height
                % 双线性插值
                pixel_value = bilinear_interpolate(img, x_img, y_img);
                corrected_img(y_corr, x_corr, :) = pixel_value;
                filled_pixels = filled_pixels + 1;
            end
        end
    end
    
    fprintf('成功填充像素: %d / %d (%.1f%%)\n', ...
            filled_pixels, width_pixels * height_pixels, ...
            filled_pixels * 100 / (width_pixels * height_pixels));
    
    % 检查输出图像是否全黑
    if filled_pixels == 0
        warning('没有像素被填充！检查单应性矩阵和坐标映射。');
    end
end

function H = compute_homography(corr_pts, img_pts)
    % H阵，矫正后/前角点坐标
    n = size(corr_pts, 1);
    A = zeros(2*n, 9);
    
    for i = 1:n
        x_corr = corr_pts(i, 1); y_corr = corr_pts(i, 2);
        x_img = img_pts(i, 1); y_img = img_pts(i, 2);
        
        A(2*i-1, :) = [x_corr, y_corr, 1, 0, 0, 0, -x_img*x_corr, -x_img*y_corr, -x_img];
        A(2*i, :)   = [0, 0, 0, x_corr, y_corr, 1, -y_img*x_corr, -y_img*y_corr, -y_img];
    end
    
    % 使用SVD求解
    [~, ~, V] = svd(A);
    H = reshape(V(:,9), 3, 3)';
    H = H / H(3,3); % 归一化
end

function pixel = bilinear_interpolate(img, x, y)
    % 双线性插值，处理边界情况
    
    x1 = floor(x);
    y1 = floor(y);
    x2 = x1 + 1;
    y2 = y1 + 1;
    
    % 边界处理
    if x1 < 1, x1 = 1; end
    if y1 < 1, y1 = 1; end
    if x2 > size(img, 2), x2 = size(img, 2); end
    if y2 > size(img, 1), y2 = size(img, 1); end
    
    % 确保坐标有效
    x1 = max(1, min(x1, size(img,2)));
    x2 = max(1, min(x2, size(img,2)));
    y1 = max(1, min(y1, size(img,1)));
    y2 = max(1, min(y2, size(img,1)));
    
    % 权重计算
    wx = x - x1;
    wy = y - y1;
    
    % 确保权重在[0,1]范围内
    wx = max(0, min(1, wx));
    wy = max(0, min(1, wy));
    
    % 插值计算
    pixel = zeros(1, 3);
    for channel = 1:3
        % 获取四个角点的值
        f11 = double(img(y1, x1, channel));
        f21 = double(img(y1, x2, channel));
        f12 = double(img(y2, x1, channel));
        f22 = double(img(y2, x2, channel));
        
        % 双线性插值公式
        top = (1 - wx) * f11 + wx * f21;
        bottom = (1 - wx) * f12 + wx * f22;
        pixel(channel) = (1 - wy) * top + wy * bottom;
    end
    
    pixel = uint8(pixel);
end