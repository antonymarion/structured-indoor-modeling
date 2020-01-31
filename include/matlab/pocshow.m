% POCSHOW - show point cloud
% POCSHOW(point, h)
% POCSHOW(point, h, MarkerSize)
% point - [X,Y,Z]
% h - figure handler
function pocshow(point, h, varargin)
nVarargs = length(varargin);
if h > 0
figure(h)
else
figure('WindowStyle','docked','Name','PointCloud') 
end
if nVarargs == 0
plot3(point(:,1), point(:,2), point(:,3), '.', 'MarkerSize', 1), axis equal
else
plot3(point(:,1), point(:,2), point(:,3), '.', 'MarkerSize', varargin{1}), axis equal  
end
end