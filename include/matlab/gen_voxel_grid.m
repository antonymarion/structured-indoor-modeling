function voxel_param = gen_voxel_grid(point_cloud, leaf_size)
% voxel_param = GENVOXELGRID(point_cloud, leaf_size)
% Generate the bounding voxel grid of the entire pointcloud

xmin = 1.0e12;
xmax = -1.0e12;
ymin = 1.0e12;
ymax = -1.0e12;
zmin = 1.0e12;
zmax = -1.0e12;

xvec = point_cloud(:,1);
xmin_ = min(xvec);
xmax_ = max(xvec);
xmin = min(xmin,xmin_);
xmax = max(xmax,xmax_);

yvec = point_cloud(:,2);
ymin_ = min(yvec);
ymax_ = max(yvec);
ymin = min(ymin,ymin_);
ymax = max(ymax,ymax_);

zvec = point_cloud(:,3);
zmin_ = min(zvec);
zmax_ = max(zvec);
zmin = min(zmin,zmin_);
zmax = max(zmax,zmax_);

% bounding box parameters
b0 = [xmin-0.05*(xmax-xmin), ymin-0.05*(ymax-ymin), zmin-0.05*(zmax-zmin)]';
gNumX = ceil(1.10*(xmax-xmin)/(leaf_size(1)));
gNumY = ceil(1.10*(ymax-ymin)/(leaf_size(2)));
gNumZ = ceil(1.10*(zmax-zmin)/(leaf_size(3)));
bSize = zeros(3,1);
bSize(1) = gNumX*leaf_size(1);
bSize(2) = gNumY*leaf_size(2);
bSize(3) = gNumZ*leaf_size(3);
gNum = [gNumX gNumY gNumZ]';
voxel_param  = struct('gNum', gNum, 'b0', b0, 'bSize', bSize);
end