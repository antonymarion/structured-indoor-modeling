%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Evidences, voxel_param] = compute_evidences(POINT, CAM, CAMLIST, core_point_idx, options, varargin)
% COMPUTE_EVIDENCES - compute evidences from the point cloud and camera information
% COMPUTE_EVIDENCES(POINT, CAM, CAMLIST, CORE_POINT_IDX, OPTIONS, PROPERTIES)
% PROPERTIES 'Property1', 'Value1', 'Property2', 'Value2',...): 
% 'PNT_ALL' : Compute Point Evidence All 
%           1 : On
%           0 : OFF
% 'PNT_WALL' : Compute Point Evidence WAll 
%           1 : On
%           0 : OFF
% 'FS_DENSE' : Compute Freespace Evidence Dense 
%           1 : On
%           0 : OFF
% 'FS_SPARSE' : Compute Freespace Evidence Sparse 
%           1 : On
%           0 : OFF
% 'NML' : Compute normal
%           1 : On
%           0 : OFF
% Theoretically, you can use different voxel grid for different structural element,
% however, in this implementation we define the unique voxel grid for the
% simplicity


Properties = struct('PNT_ALL',1, 'PNT_WALL', 1, 'FS_DENSE', 1, 'FS_SPARSE', 1, 'NML', 1); % default setting

if nargin > 3
    if mod(nargin-3, 2) ~= 0, error('Invalid options'); end
    for i = 1 : 2 : length(varargin)
            Property = varargin{i};								% read one line from file
            Value = varargin{i+1};

            switch Property
                
                case 'PNT_ALL'
                    if Value == 1
                        Properties.PNT_ALL = 1;
                    elseif Value == 0
                        Properties.PNT_ALL = 0;            
                    else
                        error(' Value of ''PNT_ALL'' must be 1 or 0 ');     
                    end
                    
                case 'PNT_WALL'
                    if Value == 1
                        Properties.PNT_WALL = 1;
                    elseif Value == 0
                        Properties.PNT_WALL = 0;            
                    else
                        error(' Value of ''PNT_WALL'' must be 1 or 0 ');     
                    end
                
               case 'FS_DENSE'
                    if Value == 1
                        Properties.FS_DENSE = 1;
                    elseif Value == 0
                        Properties.FS_DENSE = 0;            
                    else
                        error(' Value of ''FS_DENSE'' must be 1 or 0 ');     
                    end  
               case 'FS_SPARSE'
                    if Value == 1
                        Properties.FS_SPARSE = 1;
                    elseif Value == 0
                        Properties.FS_SPARSE = 0;            
                    else
                        error(' Value of ''FS_SPARSE'' must be 1 or 0 ');     
                    end 
                    
               case 'NML'
                    if Value == 1
                        Properties.NML = 1;
                    elseif Value == 0
                        Properties.NML = 0;            
                    else
                        error(' Value of ''NML'' must be 1 or 0 ');     
                    end 
                    
                case 'VIS'
                    if Value == 1
                        Properties.VIS = 1;
                    elseif Value == 0
                        Properties.VIS = 0;            
                    else
                        error(' Value of ''VIS'' must be 1 or 0 ');     
                    end 
                    
                    
                    
                otherwise
                    error(' ''%s'': Invalid property', Property);                            
            end          
    end    
end

display('[EVIDENCE COMPUTATION]')
% compute the evidences for entier points
P = cell2mat(POINT);
C = cell2mat(CAM);

wall_indeces = find(abs(P(:,11))<0.01);

if ~isempty(core_point_idx)

Pwall = P(intersect(core_point_idx , wall_indeces),:);
P = P(core_point_idx ,:);
C = C(core_point_idx );
else
Pwall = P(wall_indeces,:);
end

XYZwall = Pwall(:,3:5);
Nwall = Pwall(:,9:11);
XYZ = P(:,3:5);

point_evidence = [];
point_evidence_wall = [];
point_evidence3d = [];
freespace_evidence = [];
freespace_evidence3d = [];
freespace_evidence_sparse = [];

nx = [];
ny = [];

%% filter the original point could in the grid
display(sprintf('#init points = %d', size(XYZ,1)));
voxel_param = gen_voxel_grid(XYZ, options.leafsize);
filtered_index = mex_vgFilter(XYZ', voxel_param.gNum, voxel_param.bSize, voxel_param.b0);
display(sprintf('#filtered points = %d', size(filtered_index,1)));
display(sprintf('Voxel Grid Size = [%d(x) %d(y) %d(z)]', voxel_param.gNum(1), voxel_param.gNum(2), voxel_param.gNum(3)));
cam_filtered = C(filtered_index);
cloud_filtered = P(filtered_index,:);

Z = floor(voxel_param.gNum(3)*(CAMLIST(:,3)-voxel_param.b0(3))/voxel_param.bSize(3))+1;

if Properties.PNT_ALL == 1
display('...Computing the point evidence (all)')
point_evidence = mex_genPointEvidence(XYZ', voxel_param.gNum, voxel_param.bSize, voxel_param.b0);
point_evidence3d = reshape(point_evidence, voxel_param.gNum(2), voxel_param.gNum(1), voxel_param.gNum(3));
if options.mor3d > 0
    point_evidence = imclose(point_evidence3d>0, strel('disk', options.mor3d));
else
    point_evidence = point_evidence3d>0;
end
point_evidence = reshape(point_evidence, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
point_evidence = sum((point_evidence>0)');
point_evidence = reshape(point_evidence, voxel_param.gNum(2), voxel_param.gNum(1));
point_evidence = normalize_evidence(point_evidence, options.alpha);
point_evidence3d = normalize_evidence(point_evidence3d, options.alpha);
end

if Properties.PNT_WALL == 1
display('...Computing the point evidence (wall)')
[point_evidence_wall nx ny nz] = mex_genPointNormalEvidence(XYZwall', Nwall', voxel_param.gNum, options.leafsize, voxel_param.b0);
point_evidence_wall = reshape(point_evidence_wall, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
for z = 1 : voxel_param.gNum(3)
point_evidence_wall(:,z) = point_evidence_wall(:,z)*(1-z/voxel_param.gNum(3));
end
point_evidence_wall = sum((point_evidence_wall)');
point_evidence_wall = reshape(point_evidence_wall, voxel_param.gNum(2), voxel_param.gNum(1));
point_evidence_wall = normalize_evidence(point_evidence_wall, options.alpha);
end


XYZ_filtered = cloud_filtered(:,3:5);

if Properties.FS_DENSE == 1
display('...Computing the free-space evidence (dense)')
freespace_evidence = mex_genFreeSpaceEvidence(XYZ_filtered', cam_filtered, CAMLIST', voxel_param.gNum, options.leafsize, voxel_param.b0, 1.0);
freespace_evidence3d = reshape(freespace_evidence, voxel_param.gNum(2), voxel_param.gNum(1), voxel_param.gNum(3));
if options.mor3d > 0
    freespace_evidence = imclose(freespace_evidence3d>0, strel('disk', options.mor3d));
else
    freespace_evidence = freespace_evidence3d>0;
end
freespace_evidence = reshape(freespace_evidence, voxel_param.gNum(2)*voxel_param.gNum(1), voxel_param.gNum(3));
freespace_evidence = sum((freespace_evidence)');
freespace_evidence = reshape(freespace_evidence, voxel_param.gNum(2), voxel_param.gNum(1));
freespace_evidence = normalize_evidence(freespace_evidence, options.alpha);
freespace_evidence3d = normalize_evidence(freespace_evidence3d, options.alpha);
end

if Properties.FS_SPARSE == 1
mask = freespace_evidence > 0;
mask = imerode(mask, strel('disk',5));
freespace_evidence_sparse = freespace_evidence;
freespace_evidence_sparse(mask==0) = 0;    
% display('...Computing the free-space evidence (sparse)')
% freespace_evidence_sparse = mex_genFreeSpaceEvidence(XYZ_filtered', cam_filtered, CAMLIST', voxel_param.gNum, options.leafsize, voxel_param.b0, 0.95);
% freespace_evidence_sparse = reshape(freespace_evidence_sparse, voxel_param.gNum(2)*voxel_param.gNum(1), voxel_param.gNum(3));
% freespace_evidence_sparse = sum((freespace_evidence_sparse>0)');
% freespace_evidence_sparse = reshape(freespace_evidence_sparse, voxel_param.gNum(2), voxel_param.gNum(1));
% freespace_evidence_sparse = normalize_evidence(freespace_evidence_sparse, options.alpha);
end

if Properties.NML == 1
display('...Computing the normal map')
nx = reshape(nx, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
nx = sum((nx)');
nx = reshape(nx, voxel_param.gNum(2), voxel_param.gNum(1));
ny = reshape(ny, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
ny = sum((ny)');
ny = reshape(ny, voxel_param.gNum(2), voxel_param.gNum(1));
nz = reshape(nz, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
nz = sum((nz)');
nz = reshape(nz, voxel_param.gNum(2), voxel_param.gNum(1));
nm = sqrt(nx.*nx + ny.*ny + nz.*nz);
nx = nx./nm;
ny = ny./nm;
nx(isnan(nx)) = 0;
ny(isnan(ny)) = 0;
end

% Evidences = struct('point_evidence', point_evidence, 'point_evidence_wall', point_evidence_wall, 'point_evidence3d', point_evidence3d, 'point_evidence_wall3d', point_evidence_wall3d, 'freespace_evidence', freespace_evidence, 'freespace_evidence3d', freespace_evidence3d, 'freespace_evidence_sparse', freespace_evidence_sparse, 'featureVis', featureVis_dense, 'nx', nx, 'ny', ny);
Evidences = struct('point_evidence', point_evidence, 'point_evidence_wall', point_evidence_wall, 'point_evidence3d', point_evidence3d, 'freespace_evidence', freespace_evidence, 'freespace_evidence3d', freespace_evidence3d, 'freespace_evidence_sparse', freespace_evidence_sparse, 'nx', nx, 'ny', ny);
end
