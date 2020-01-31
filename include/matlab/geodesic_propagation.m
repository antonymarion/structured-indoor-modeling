function propagated = geodesic_propagation(src, BW)
%% propagate room labels to others
BWcopy = BW;
[h w] = size(BW);
BW(find(BW==0)) = rand(size(find(BW==0)));
[dIf dIb idf idb] = calc_dIandId_gray(BW(:), h, w);
M = zeros(h*w,1);
L = zeros(h*w,1);
idzero = find(src~=0);

L0 = [1:h*w]';
L0 = reshape(L0,h,w);
[X, Y] = meshgrid([1:w],[1:h]);
id = (X(:)-1)*h + Y(:);
M0 = 1.0e12*ones(h*w,1);
M0(idzero) = 0;
M0 = reshape(M0,h,w);
% geodesic distance transform
r = 1;
lambda = 1.0e12;
[M, L] = mex_gdt(M0, L0, dIf, dIb, r, lambda, 4);
M = reshape(M,h,w);
L = reshape(L,h,w);
propagated = src(L);
propagated(BWcopy==0) = 0;
propagated = reshape(propagated, h,w);
end

