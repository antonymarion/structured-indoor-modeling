%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [floorHeight, ceilHeight] = floor_ceil_reconstruction(POINT_room, CAM_room, CAMLIST)
% [floorHeight, ceilHeight] = FLOOR_CEIL_RECONSTRUCTION(POINT_room, CAM_room, CAMLIST)
% RANSAC-based plane fitting for floor and ceiling extraction

%% (b) floor segmentation
floorHeight = zeros(size(POINT_room,1),1);
ceilHeight = zeros(size(POINT_room,1),1);

%% Extract points whose normals are directed at [0 0 1] or [0 0 -1]
PFLOORCEIL = cell(1,1);
CFLOORCEIL = cell(1,1);
for i = 1 : length(POINT_room)
    PFLOORCEIL{i,1} = POINT_room{i}(find(abs(POINT_room{i}(:,11)) > 0.9),:);
    CFLOORCEIL{i,1} = CAM_room{i}(find(abs(POINT_room{i}(:,11)) > 0.9),1);
end


floorParam = cell(1,1);
ceilParam = cell(1,1);
for i = 1 : length(POINT_room)
FP = room_floor_extraction(PFLOORCEIL{i}, CFLOORCEIL{i}, CAMLIST);

if ~isempty(FP)
floorParam{i,1} = FP;
[XX,YY] = meshgrid([min(PFLOORCEIL{i}(:,3)):(max(PFLOORCEIL{i}(:,3))-min(PFLOORCEIL{i}(:,3)))/100:max(PFLOORCEIL{i}(:,3))],[min(PFLOORCEIL{i}(:,4)):(max(PFLOORCEIL{i}(:,4))-min(PFLOORCEIL{i}(:,4)))/100:max(PFLOORCEIL{i}(:,4))]);
ZZ = -(floorParam{i}(1)*XX + floorParam{i}(2)*YY + floorParam{i}(4))/floorParam{i}(3);

figure('WindowStyle', 'docked', 'Name', sprintf('RANSAC:Floor(%d)', i)),
plot3(PFLOORCEIL{i}(:,3),PFLOORCEIL{i}(:,4),PFLOORCEIL{i}(:,5),'.', 'MarkerSize',1)
hold on
plot3(XX,YY,ZZ)
axis equal
hold off
title(mean(ZZ(:)))
floor_height = mean(ZZ(:));
floorHeight(i) = floor_height;
else
floor_hegith(i) = median(floorHeight(:));    
end

CP = room_ceiling_extraction(PFLOORCEIL{i}, CFLOORCEIL{i}, CAMLIST);

if ~isempty(CP)
ceilParam{i,1} = CP; 
[XX,YY] = meshgrid([min(PFLOORCEIL{i}(:,3)):(max(PFLOORCEIL{i}(:,3))-min(PFLOORCEIL{i}(:,3)))/100:max(PFLOORCEIL{i}(:,3))],[min(PFLOORCEIL{i}(:,4)):(max(PFLOORCEIL{i}(:,4))-min(PFLOORCEIL{i}(:,4)))/100:max(PFLOORCEIL{i}(:,4))]);
ZZ = -(ceilParam{i}(1)*XX + ceilParam{i}(2)*YY + ceilParam{i}(4))/ceilParam{i}(3);

figure('WindowStyle', 'docked', 'Name', sprintf('RANSAC:Ceil(%d)', i)),
plot3(PFLOORCEIL{i}(:,3),PFLOORCEIL{i}(:,4),PFLOORCEIL{i}(:,5),'.', 'MarkerSize',1)
hold on
plot3(XX,YY,ZZ)
axis equal
hold off
title(mean(ZZ(:)))
ceilHeight(i) = mean(ZZ(:));
else
ceilHeight(i) = median(ceilHeight(:));
end
end


end