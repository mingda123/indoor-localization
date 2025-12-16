function [anchors, epochs] = UWB_fileread(filename)
% 解析UWB定位输出数据并按历元输出
% 输入:
%   filename - 数据文件名
% 输出:
%   anchors - 基站坐标信息
%   range_measurements - 距离测量数据
%   position_estimates - 定位坐标估计
%   smoothed_positions - 平滑后的坐标

    % 初始化数据结构
    anchors = [];
    epochs = struct('timestamp', {}, 'RR', {}, 'LE', {});
    
    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end
    
    last_epoch = [];
    epoch_count = 0;
    timestamp_tolerance = 50;
    anchor_count = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line) || ~startsWith(line, 'T:')
            continue;
        end
        
        % 分割行数据
        parts = strsplit(line, ':');
        
        % 处理基站坐标信息 (AP行)
        if length(parts) >= 7 && strcmp(parts{3}, 'AP')
            anchor_count = anchor_count + 1;
            anchors(anchor_count).timestamp = str2double(parts{2});
            anchors(anchor_count).id = str2double(parts{4});
            anchors(anchor_count).x = str2double(parts{5});
            anchors(anchor_count).y = str2double(parts{6});
            anchors(anchor_count).z = str2double(parts{7});
            
        % 处理距离测量数据 (RR行)
        elseif length(parts) >= 8 && strcmp(parts{3}, 'RR')
            timestamp = str2double(parts{2});
            anchor_id = int32(str2double(parts{5}));

            % 新历元，保存上一历元
            if ~isempty(last_epoch) && abs(timestamp - last_epoch.timestamp) > timestamp_tolerance
                epoch_count = epoch_count + 1;
                epochs(epoch_count) = last_epoch;
                last_epoch = [];
            end

            % 初始化当前历元
            if isempty(last_epoch)
                last_epoch.timestamp = timestamp;
                last_epoch.RR = struct('raw_range', {}, 'corrected_range', {});
                last_epoch.LE = [];
            end
            
            % 添加距离测量
            range_data.raw_range = str2double(parts{6});
            range_data.corrected_range = str2double(parts{7});
            last_epoch.RR(anchor_id+1) = range_data;
            
        % 处理定位坐标估计 (LE行)
        elseif length(parts) >= 9 && strcmp(parts{3}, 'LE')
            timestamp = str2double(parts{2});
            
            % 如果时间戳匹配当前历元，添加定位结果
            if ~isempty(last_epoch) && abs(timestamp - last_epoch.timestamp) < timestamp_tolerance
                coord_str = parts{7};
                coord_str = strrep(coord_str, '[', '');
                coord_str = strrep(coord_str, ']', '');
                coords = strsplit(coord_str, ',');
                
                last_epoch.LE.x = str2double(coords{1});
                last_epoch.LE.y = str2double(coords{2});
                last_epoch.LE.z = str2double(coords{3});
            end
            
        end
    end
    
    fclose(fid);
    
    % 显示解析结果统计
    fprintf('数据解析完成 %s:\n',filename);
    fprintf('  基站数量: %d\n', anchor_count);
    fprintf('  历元数: %d\n', epoch_count);
end