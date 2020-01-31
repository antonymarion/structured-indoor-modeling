function [options, core_point_idx] = point_noise_removal(POINT, options)
% [options, core_point_idx] = POINT_NOISE_REMOVAL(POINT, options)
% Remove noises in input point cloud

display('...computing reliable point indeces');

if options.noiseremoval == 0
% Use only intensity information in the input pointcloud for removing
% unreliable points

intensity_thresh = zeros(length(POINT),1);
core_point_idx = cell(1,1);
count = 0;
for i = 1 : length(POINT)
    P = POINT{i};
    intensity_thresh(i) = mean(P(:,12)) - options.std*std(P(:,12)); % default 2   
    core_point_idx{i,1} = find(P(:,12)>=intensity_thresh(i)) + count;
    count = count + size(P,1);
end
options.intensity_thresh = intensity_thresh;
core_point_idx = cell2mat(core_point_idx);
clear P intensity_thresh count;

else
% Use geometric distribution of the input pointcloud for removing
% unreliable points
    
core_point_idx = cell(1,1);
count = 0;
for i = 1 : length(POINT)
fprintf('*');
P = POINT{i};
XYZ = P(:,3:5);
intensity_thresh = mean(P(:,12)) - 1.5*std(P(:,12)); % default 2   
index_intensity = find(P(:,12)>=intensity_thresh);

voxel_param = gen_voxel_grid(XYZ, options.leafsize*4);
point_evidence = mex_genPointEvidence(XYZ', voxel_param.gNum, voxel_param.bSize, voxel_param.b0);
point_evidence = reshape(point_evidence, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
point_evidence = sum((point_evidence>0)');
point_evidence = reshape(point_evidence, voxel_param.gNum(2), voxel_param.gNum(1));
mask = point_evidence > 0;
se = strel('disk',1);
mask = imerode(mask, se);
mask = imdilate(mask, se);
CC = bwconncomp(mask,8);
numPix = zeros(length(CC.PixelIdxList),1);
for k = 1 : length(CC.PixelIdxList)
    numPix(k) = size(CC.PixelIdxList{k},1);
end
[~, maxidx] = max(numPix);
mask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
mask(CC.PixelIdxList{maxidx}) = 1;
B = bwboundaries(mask,8);
boundary_img = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
boundary_img((B{1}(:,2)-1)*voxel_param.gNum(2)+B{1}(:,1)) = 1;
mask = imfill(boundary_img,'holes');

% project 3d to 2d
X = floor(voxel_param.gNum(1)*(XYZ(:,1)-voxel_param.b0(1))/voxel_param.bSize(1))+1;
Y = floor(voxel_param.gNum(2)*(XYZ(:,2)-voxel_param.b0(2))/voxel_param.bSize(2))+1; 

idx = (X-1)*voxel_param.gNum(2) + Y;
core_point_idx{i,1} = intersect(find(mask(idx)==1),index_intensity) + count;
count = count + size(XYZ,1);
end
fprintf('\n');
core_point_idx = cell2mat(core_point_idx);
end
close all
end
