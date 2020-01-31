function [dIf dIb idf, idb] = calcdIandId_gray(img, h, w)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述

dIf = zeros(h*w,5);       
dIb = zeros(h*w,5);

% forward for Vf(x,y)
% Vf(0,0)
dIf(:,1) = 0;      

% Vf(-1,-1)
[xg1 yg1] = meshgrid([1:w-1],[1:h-1]); % neighbor
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([2:w],[2:h]); % target
id2 = (xg2-1)*h + yg2;
dIf(id2,2) = (img(id1,1)-img(id2,1)).^2;

% Vf(0,-1)
[xg1 yg1] = meshgrid([1:w],[1:h-1]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([1:w],[2:h]);
id2 = (xg2-1)*h + yg2;
dIf(id2,3) = (img(id1,1)-img(id2,1)).^2;

% Vf(1,-1)
[xg1 yg1] = meshgrid([2:w],[1:h-1]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([1:w-1],[2:h]);
id2 = (xg2-1)*h + yg2;
dIf(id2,4) = (img(id1,1)-img(id2,1)).^2;


% Vf(-1,0)
[xg1 yg1] = meshgrid([1:w-1],[1:h]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([2:w],[1:h]);
id2 = (xg2-1)*h + yg2;
dIf(id2,5) = (img(id1,1)-img(id2,1)).^2;
dIf = sqrt(dIf);


% backward for Vb(x,y)
% Vb(0,0)
dIb(:,1) = 0;      

% Vb(1,1)
[xg1 yg1] = meshgrid([2:w],[2:h]); % neighbor
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([1:w-1],[1:h-1]); % target
id2 = (xg2-1)*h + yg2;
dIb(id2,2) = (img(id1,1)-img(id2,1)).^2;

% Vf(0,1)
[xg1 yg1] = meshgrid([1:w],[2:h]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([1:w],[1:h-1]);
id2 = (xg2-1)*h + yg2;
dIb(id2,3) = (img(id1,1)-img(id2,1)).^2;


% Vf(-1,1)
[xg1 yg1] = meshgrid([1:w-1],[2:h]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([2:w],[1:h-1]);
id2 = (xg2-1)*h + yg2;
dIb(id2,4) = (img(id1,1)-img(id2,1)).^2;


% Vf(1,0)
[xg1 yg1] = meshgrid([2:w],[1:h]);
id1 = (xg1-1)*h + yg1;
[xg2 yg2] = meshgrid([1:w-1],[1:h]);
id2 = (xg2-1)*h + yg2;
dIb(id2,5) = (img(id1,1)-img(id2,1)).^2;
dIb = sqrt(dIb);

Vfx = [0,-1,0,1,-1];
Vfy = [0,-1,-1,-1,0];

% 
[X, Y] = meshgrid([1:w],[1:h]);
out = ones(h*w,5);
X = X(:);
Y = Y(:);
X_ = repmat(X,1,5)+repmat(Vfx,h*w,1);
Y_ = repmat(Y,1,5)+repmat(Vfy,h*w,1);
out(X_ < 1 | X_ > w) = 0;
out(Y_ < 1 | Y_ > h) = 0;
idf = (X_-1)*h+Y_;
idf(out==0) = -1;

Vbx = [0,1,0,-1,1];
Vby = [0,1,1,1,0];
out = ones(h*w,5);
X_ = repmat(X,1,5)+repmat(Vbx,h*w,1);
Y_ = repmat(Y,1,5)+repmat(Vby,h*w,1);
out(X_ < 1 | X_ > w) = 0;
out(Y_ < 1 | Y_ > h) = 0;
idb = (X_-1)*h+Y_;
idb(out==0) = -1;

end

