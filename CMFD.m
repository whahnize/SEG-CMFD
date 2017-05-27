
%% Initialization and Adding library
DEBUG = true;

if DEBUG == true
    path = 'data/MICC_F600/red_tower.png';
else
    clc; clear; 
    [file, full_path] = uigetfile('*', 'Pick  a picture to be inspected copy-move forgery ');
    path = strcat(full_path, file);
end

% Add vlfeat library: vl_* 
addpath(genpath('vlfeat-0.9.20-bin'));

%% Read an image and transform into single type
copy_img = imread(path);
Is = im2single(copy_img);
[c_r, c_c, c_channel] = size(Is);
%% Parameter settings for SLIC-based segmentation (parameters are chosen for MICC_F600 dataset)

split_path = strsplit(path,'/');
set = split_path(2);

if strcmpi(set,'MICC_F600')
    % parameters for vl_slic
    r_size = 50;
    reg = 0.8;
    img_size = c_r * c_c;

    if img_size > 3000*2000
        r_size = 200;
    elseif img_size > 2000*1000
        r_size = 150;
    elseif img_size > 1000*600
        r_size = 100;
    else
        r_size = 50;
    end
    SEGMENTS = vl_slic(Is, r_size, reg)+1;
elseif strcmpi(set, 'Christlein')
    % parameters for vl_quickseg
    ratio = 0.7;
    kernel = 1;
    % TO DO: figure out what parameter 'MAXDIST' means and what an adequate value is  
    maxdist = 0;
    
    SEGMENTS = vl_quickseg(Is, 0.7,1,maxdist)+1;
else
    disp('dataset cannot be identified');
    return
end

%% Visualize segmentation result
% VisSegmentation(copy_img,SEGMENTS);

%% First Matching (Robust matching)
%% 1. Keypoint Extraction (SIFT)
% Compute gray scale image
Ig = rgb2gray(copy_img);

% Compute SIFT features
[f,d] = vl_sift(single(Ig));

% Patch - keypoints map (patch_num: [descriptor of keypoint1,descriptor of keypoint2, ...])
patchToKeypoints = containers.Map('KeyType','int32','ValueType','any');

% Initialize container
maxPatch = double(max(SEGMENTS(:)));
for k = 1:maxPatch
    if ~(isKey(patchToKeypoints, k))
        patchToKeypoints(k) = [];
    end
end

% Normalize descriptors and assign them to each patch
[feature_elements, num_of_features] = size(f);
norm_d = zeros(size(d));
for n_f = 1:num_of_features
    feature_x = f(1,n_f);
    feature_y = f(2,n_f);
    feature_d = double(d(:,n_f));
    % first norm
    feature_d = feature_d / norm(feature_d);
    % trancate using 0.2 as a threshold 
    feature_d(feature_d>0.2) = 0.2;
    % second norm
    feature_d = feature_d / norm(feature_d);
    norm_d(:,n_f) = feature_d;
    cur_seg = SEGMENTS(int32(feature_y),int32(feature_x));
    patchToKeypoints(cur_seg) = [patchToKeypoints(cur_seg); transpose(feature_d)];
end

% Patch matching using k-d tree
kdtree = vl_kdtreebuild(norm_d);
dis_threshold = 0.4;
for p1 = 1:1
    keypoints_p1 = patchToKeypoints(p1); % [descriptor 1; descriptor 2; ...]
    num_of_keypoints_p1 = size(keypoints_p1, 1);
    point1_set = transpose(keypoints_p1);
    num_matched_keypoints = zeros(maxPatch,1);
    for keypoint = 1:num_of_keypoints_p1
        % k-nearest neighbers 11= 1 + 10
        [index, distance] = vl_kdtreequery(kdtree, norm_d, point1_set(:, keypoint), 'NumNeighbors', 11); 
        thres_matched_points = floor(f(1:2,index(distance<dis_threshold)));
        x_coords = thres_matched_points(1,:);
        y_coords = thres_matched_points(2,:);
        matched_patches = diag(SEGMENTS(y_coords,x_coords));
        % ISSUE: should we include its own patch or exclude it when setting
        % threshold?
        
        % Add 1 to the patches whose descriptors are matched by query
        % descriptor including own patch
        % num_matched_keypoints(matched_patches) = num_matched_keypoints(matched_patches) + 1;
        
        % Excluded matching
        num_matched_keypoints(matched_patches~=p1) = num_matched_keypoints(matched_patches~=p1) + 1;
    end
    threshold = 10 * sum(num_matched_keypoints) / maxPatch; 
    suspicious_patches = find(num_matched_keypoints>threshold & num_matched_keypoints~=p1);
    
    if ~isempty(suspicious_patches)
        X = sprintf(' index of source patch: %d\n num of suspicious patches: %d',p1,sum(suspicious_patches));
        disp(X);
    end
    
end
% matched_list = [];
% for p1 = 0:maxPatch-1
%     keypoints_p1 = patchToKeypoints(p1); % [descriptor 1; descriptor 2; ...]
%     num_of_keypoints_p1 = size(keypoints_p1, 1);
%     point1_set = transpose(keypoints_p1);
%     for p2 = p1+1:maxPatch
%         num_matched_keypoints = 0;
%         keypoints_p2 = patchToKeypoints(p2); % [descriptor 1; descriptor 2; ...]
%         if isempty(keypoints_p2) 
%             continue 
%         end
%         num_of_keypoints_p2 = size(keypoints_p2, 1);
%         point2_set = transpose(keypoints_p2);
%         kdtree = vl_kdtreebuild(point2_set);
%         [index,distance] = vl_kdtreequery(kdtree, point2_set, point1_set, 'NumNeighbors', 10);
%         
%         for k = 1:num_of_keypoints_p1
%             for cand = 1:10
%                 if distance(cand, k) < 0.16
%                     num_matched_keypoints = num_matched_keypoints + 1;
%                 end

%             end
%         end
%         if num_matched_keypoints > threshold
%             display(num_matched_keypoints);
%         end
%     end
% end

%% Second Matching (Iteration)
%% 
