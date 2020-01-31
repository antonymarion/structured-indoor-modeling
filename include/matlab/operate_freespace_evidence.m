function [BW2_resized, BW2] = operateFreespaceEvidence(BW, scale, mor_r, mor_e, min_length)



BW = BW>0;

BW = imclose(BW, strel('disk', mor_r));

invmask = 1-BW;

BW2 = BW;
CC = bwconncomp(invmask);
for i = 1 : CC.NumObjects
   region_indeces = CC.PixelIdxList{i};
   if length(region_indeces)<mor_e
      BW2(region_indeces) = 1; 
   end
end

gNum = size(BW2);
theta = [0,pi/2,pi,3*pi/2]; % we can also find the manhattan direction here
vec = [cos(theta);sin(theta)];

% 1st operation
x = zeros(gNum(1),gNum(2));
y = zeros(gNum(1),gNum(2));
for i = 1 : size(vec,2);
    F = mex_feature2D(gNum, BW2, vec(:,i));
    F = reshape(F,gNum(1),gNum(2));
    if(i==1 | i==3)
        x = x + F;
    else
        y = y + F;
    end
end
z = min(x,y);
d = min_length; % mobile void 2nd

count = 0;
while(1)
    count = count + 1;
BW2_resized = z > d;
% 2nd operation
x = zeros(gNum(1),gNum(2));
y = zeros(gNum(1),gNum(2));
for i = 1 : size(vec,2);
    F = mex_feature2D(gNum, BW2_resized, vec(:,i));
    F = reshape(F,gNum(1),gNum(2));
    if(i==1 | i==3)
        x = x + F;
    else
        y = y + F;
    end
end
z = min(x,y);

BW2_new = z > d;
thresh = 0.00000001*length(find(BW2_resized>0));
if sum(sum(abs(double(BW2_new)-double(BW2_resized)))) <= 0
    BW2_resized = BW2_new;
    break;
end
BW2_resized = BW2_new;
end

%% downsampling
[X Y] = meshgrid([1:scale:gNum(2)], [1:scale:gNum(1)]);
index = (X(:)-1)*gNum(1) + Y(:);
BW2_resized = imresize(BW2_resized, 1/scale);
BW2_resized = double(BW2_resized>0.95);
end

