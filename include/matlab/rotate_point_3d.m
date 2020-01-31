function [newPOINT] = rotatePoints(POINT, rotmat)   

if iscell(POINT)
numCluster = size(POINT,1);
newPOINT = POINT;
for i = 1 : numCluster
    P = POINT{i,1};
    XYZ = (rotmat*P(:,3:5)')';
    P(:,3:5) = XYZ;
    NXYZ = (rotmat*P(:,9:11)')';
    P(:,9:11) = NXYZ;
    newPOINT{i,1} = P;
end

else
    
    if size(POINT,2) == 12
    P = POINT;
    XYZ = (rotmat*P(:,3:5)')';
    P(:,3:5) = XYZ;
    NXYZ = (rotmat*P(:,9:11)')';
    P(:,9:11) = NXYZ;
    newPOINT = P;
    end
    
    if size(POINT,2) == 3
    P = POINT;
    newPOINT = (rotmat*P(:,1:3)')';
    end
    
    
end


end

