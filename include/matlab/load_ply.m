function [POINT, CAM, CAMLIST] = load_ply(datapath)
% [POINT, CAM, CAMLIST] = LOAD_PLY(datapath)
% load ply files 
% datapath : the path of input directory
% POINT: store the point information 
% CAM: store the camera index a point is captured
% CAMLIST: store the camera position

plyList = dir(sprintf('%s/*.ply', datapath));
numPly = length(plyList); 

% load point cloud and camera data
CAMLIST = zeros(numPly,3);
POINT = cell(numPly,1);
CAM = cell(numPly,1);

for i = 1 : numPly
fprintf('loading %s\n', sprintf('%s/%s', datapath, plyList(i).name));
[cam, point] = read_ply(sprintf('%s/%s', datapath, plyList(i).name));
CAMLIST(i,:) = cam;
POINT{i} = point;
CAM{i} = i*ones(size(point,1),1);
clear point;
clear cam;
end



