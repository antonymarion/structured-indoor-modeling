% mex_plane_fitting_ransac has dependency on Eigen3.2.2
% The example of the mex compilation script should be
% mex mex_plane_fitting_ransac.cpp -I'[DIR]\Eigen3.2.2'
function maxParam = room_ceiling_extraction(P, C, CAMLIST)
% maxParam = ROOM_CEILING_EXTRACTION(P, C, CAMLIST)
camlist = unique(C);
Pnew = [];
for k = 1 : length(camlist)
camIndex = find(C==camlist(k));
Psub = P(camIndex,3:5);
index = find(Psub(:,3) >= CAMLIST(camlist(k),3));
Psub = Psub(index,:);
Pnew = [Pnew;Psub];
end

n0 = [0,0,1];
t = 1.0e-8;
if size(Pnew,1) > 100
[maxParam, inliers, iter] = mex_plane_fitting_ransac(Pnew(:,1), Pnew(:,2), Pnew(:,3), n0, 0.01*pi, t);
else
    maxParam = [];
end
end

