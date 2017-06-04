
%% Initialization and Adding library
clear all; clc;
DEBUG = true;

if DEBUG == true
    path = 'data/MICC_F600/horses.png';
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
seg_im=VisSegmentation(copy_img,SEGMENTS);
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
dis_threshold = realmax;
patch_pair = [];
patch_keypoints = [];
for p1 = 1:maxPatch
    keypoints_p1 = patchToKeypoints(p1); % [descriptor 1; descriptor 2; ...]
    num_of_keypoints_p1 = size(keypoints_p1, 1);
    point1_set = transpose(keypoints_p1);
    num_matched_keypoints = zeros(maxPatch,1);
    matched_keypoints = cell(1,maxPatch);
    for keypoint = 1:num_of_keypoints_p1
        % k-nearest neighbers 11= 1 + 10
        [index, distance] = vl_kdtreequery(kdtree, norm_d, point1_set(:, keypoint), 'NumNeighbors', 11); 
        thres_matched_points = f(1:2,index(distance<dis_threshold));
        x_coords = floor(thres_matched_points(1,:));
        y_coords = floor(thres_matched_points(2,:));
        matched_patches = diag(SEGMENTS(y_coords,x_coords));
        % ISSUE: should we include its own patch or exclude it when setting
        % threshold?
        
        % Included matching
        % Add 1 to the patches whose descriptors are matched by query
        % descriptor including own patch
        num_matched_keypoints(matched_patches) = num_matched_keypoints(matched_patches) + 1;
        
        % Store matched keypoints for each patches
        for patch_index = 1:size(matched_patches,1)
            matched_keypoints{1,matched_patches(patch_index)} = [matched_keypoints{1,matched_patches(patch_index)}; thres_matched_points(:,1).' thres_matched_points(:,patch_index).'];
        end
        
        % Excluded matching
        % num_matched_keypoints(matched_patches~=p1) = num_matched_keypoints(matched_patches~=p1) + 1;
    end
    match_threshold = 10 * sum(num_matched_keypoints) / maxPatch; 
    suspicious_patches = find(num_matched_keypoints>match_threshold);
    % Exclude src patch and replicated pair
    suspicious_patches = suspicious_patches(suspicious_patches>p1);
    
    if ~isempty(suspicious_patches)
        X = sprintf(' index of source patch: %d\n number of suspicious patches: %d',p1,size(suspicious_patches,1));
        disp(X);
        % p1 pair [src_patch suspicious_patch threshold num_of_p1_feature num_of_supicious_feature]
        p1_pair = [ones(size(suspicious_patches,1),1)*p1 suspicious_patches ones(size(suspicious_patches,1),1)*match_threshold ones(size(suspicious_patches,1),1)*num_matched_keypoints(p1) num_matched_keypoints(suspicious_patches)];
        patch_pair = [patch_pair; p1_pair];
        patch_keypoints = [patch_keypoints matched_keypoints(1,suspicious_patches)];
    end
    
end

% Estimate RANSAC affine transformation using matched feature points between pairs
% You can get transformed image by using 
% t_s = imwarp(image,tform,'OutputView',imref2d(size(image)));

tforms = [];
for match = 1:size(patch_keypoints,2)
    % matrixes which transform source segment to suspicious segment
    t = estimateGeometricTransform(patch_keypoints{1,match}(:,3:4),patch_keypoints{1,match}(:,1:2),'affine');
    tforms = [tforms; t];
end

% visualize suspicous patch pair
num_pair = size(patch_pair,1);
src= single(zeros(size(SEGMENTS)));
des= single(zeros(size(SEGMENTS)));
for p = 1:num_pair
    src(SEGMENTS==patch_pair(p,1)) = 1;
    des(SEGMENTS==patch_pair(p,2)) = 1;
end
B = imoverlay(seg_im,src,'blue');
imshow(B);
B = imoverlay(B,des,'red');
imshow(B);

%% Second Matching (Iteration)

% EM algorithm for estimating tform more accurately 

for pair = 1:num_pair
    src_idx = patch_pair(pair,1); 
    cmf_idx = patch_pair(pair,2);
    
    % First step: find new coordnates in transformed image whose DENSE sift
    % descriptor are more similar to copying source patch than old coordinates.
    I_src = single(Ig);
    t_segments = imwarp(SEGMENTS,invert(tforms(pair)),'OutputView',imref2d(size(SEGMENTS)));
    I_hat = imwarp(I_src,invert(tforms(pair)),'OutputView',imref2d(size(I_src)));
    
    [src_row,src_col] = find(SEGMENTS==src_idx);
    [cmf_row,cmf_col] = find(t_segments==cmf_idx);
end
%% 
