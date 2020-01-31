%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0, x1, mincost, flag] = wall_profile_analysis(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, thresh, roomid1, roomid2, wallid1, wallid2)
% [x0, x1, mincost, flag] = WALL_PROFILE_ANALYSIS(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, thresh)
% Clasify the wall connection type from the wall profiles by looking up the classification table

lambda = 0.01;

integralfsp0 = integralImage(freeSpaceProfile0);
integralfsp1 = integralImage(freeSpaceProfile1);

ib0 = integralImage(freeSpaceProfile0>0);
ib1 = integralImage(freeSpaceProfile1>0);

integralwpf0 = integralImage(wallProfile0);
integralwpf1 = integralImage(wallProfile1);

[h,w] = size(freeSpaceProfile0);

ymin = 0.5*h;
xmin = 0.5*(max_x-min_x);
margin = 0;

[x mincost, c1, c2] = mex_door_cost(integralfsp0, integralfsp1, integralwpf0, integralwpf1, freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, lambda, min_x, max_x, 0); % x is c++ coordinate -1

flag = 0;
x0 = [x(1)-1, max(x(2)-1,1)];
x1 = [x(3)+1, x(4)+1];  

if mincost < 0 & c1 < 0 & c2 < 0 & abs(x(1) - x(3)) > thresh
    
flag = -1;

[x mincost] = mex_door_cost(integralfsp0, integralfsp1, integralwpf0, integralwpf1, freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, lambda, min_x, max_x, 1);

xmin = min(x(3), x(1));
xmax = max(x(3), x(1));
[Hwp1, Fwp1] = mex_find_ceil_line(wallProfile0, freeSpaceProfile0, margin, ymin, xmin, xmax);
[Hwp2, Fwp2] = mex_find_ceil_line(wallProfile1, freeSpaceProfile1, margin, ymin, xmin, xmax);
Fwp = min(Fwp1, Fwp2); % floor height is very unreliable
[Rwp1, Lwp1] = mex_find_side_line(wallProfile0, freeSpaceProfile0, margin, xmin, Fwp1+0.1*(Hwp1-Fwp1), Hwp1-0.1*(Hwp1-Fwp1));
[Rwp2, Lwp2] = mex_find_side_line(wallProfile1, freeSpaceProfile1, margin, xmin, Fwp2+0.1*(Hwp2-Fwp2), Hwp2-0.1*(Hwp2-Fwp2));
Hwp = max(Hwp1, Hwp2);

x0 = [x(1)-1, max(x(2)-1,1)];
x1 = [x(3)+1, x(4)+1];    
xx = [x0(1)+1,x1(1)+1,x1(1)+1,x0(1)+1,x0(1)+1];
yy = [x0(2)+1,x0(2)+1,x1(2)+1,x1(2)+1,x0(2)+1];
test = [wallProfile0 wallProfile1;freeSpaceProfile0  freeSpaceProfile1];
figure('WindowStyle', 'docked', 'Name', sprintf('c1:%f c2%f x(1), x(2)', c1, c2)), imagesc(test)
axis equal
hold on
plot(xx,yy,'o-w','LineWidth',2)
plot(xx+w,yy,'o-w','LineWidth',2)
plot(xx,yy+h,'o-w','LineWidth',2)
plot(xx+w,yy+h,'o-w','LineWidth',2)
plot([1,w],[Hwp1+1,Hwp1+1],'--g','LineWidth',2)
plot([w+1,2*w],[Hwp2+1,Hwp2+1],'--g','LineWidth',2)
plot([1,w],[Fwp+1,Fwp+1],'--g','LineWidth',2)
plot([w+1,2*w],[Fwp+1,Fwp+1],'--g','LineWidth',2)

plot([Lwp1+1,Lwp1+1],[1 h],'--g','LineWidth',2)
plot([Lwp2+1+w,Lwp2+1+w],[1 h],'--g','LineWidth',2)
plot([Rwp1+1,Rwp1+1],[1 h],'--g','LineWidth',2)
plot([Rwp2+1+w,Rwp2+1+w],[1 h],'--g','LineWidth',2)

% sufficient condition (door: some marging between door and ceiling or
% there are some wings)
wing_thresh = 10;

wing1 = 3;
if Rwp1-max(x(3),x(1)) >= wing_thresh || min(x(3),x(1))-Lwp1 >= wing_thresh
    wing1 = 2;
end

if Rwp1-max(x(3),x(1)) >= wing_thresh && min(x(3),x(1))-Lwp1 >= wing_thresh
    wing1 = 1;
end

wing2 = 3;
if Rwp2-max(x(3),x(1)) >= wing_thresh || min(x(3),x(1))-Lwp2 >= wing_thresh
    wing2 = 2;
end

if Rwp2-max(x(3),x(1)) >= wing_thresh && min(x(3),x(1))-Lwp2 >= wing_thresh
    wing2 = 1;
end
    
top1 = 2;
if (Hwp1+1) - abs(max(x(4),x(2))) >= 0.05*(Hwp1-Fwp1)
    top1 = 1;
end

top2 = 2;
if (Hwp2+1) - abs(max(x(4),x(2))) >= 0.05*(Hwp2-Fwp2)
    top2 = 1;
end

conn_type = [1 2;3 4;5 6];
conn_matrix = [1 1 1 1 1 1;1 1 1 1 2 2;1 1 1 1 1 2; 1 1 1 2 2 2; 1 2 1 2 2 2; 1 2 2 2 2 2]; % connection type classification table
flag = conn_matrix(conn_type(wing1,top1), conn_type(wing2, top2));
if flag == 1 && abs(min(x(4),x(2))-(Fwp+1))<=20 && abs(x(2) - x(4)) > 2*thresh
        title(sprintf('DOOR: flag=%d (%d, %d)',flag, conn_type(wing1,top1), conn_type(wing2,top2)));
        display(sprintf('DOOR connection detected (RID1, WID1) = (%d, %d), (RID2, WID2) = (%d, %d)', roomid1, wallid1, roomid2, wallid2));
        return;
end

if (flag==1 || flag == 2) && abs(x(2) - x(4)) > 2*thresh
        title(sprintf('MERGE: flag=%d (%d, %d)',flag, conn_type(wing1,top1), conn_type(wing2,top2)));
        display(sprintf('MERGE connection detected (RID1, WID1) = (%d, %d), (RID2, WID2) = (%d, %d)', roomid1, wallid1, roomid2, wallid2));
        return;
end

title(sprintf('flag=%d',flag));
end