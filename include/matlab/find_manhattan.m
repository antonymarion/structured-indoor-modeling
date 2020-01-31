function M = find_manhattan(point, normal)
% M = FIND_MANHATTAN(point, normal)
% Find the Manhattan-World directions (dominant XYZ axis)


nt = normal(normal(:,3)<-0.9,:);
n1 = mex_ransac_manhattan_1st(nt'); % firstly get the dominant coordinate (should be z-axis)
normal_sub = normal(abs(normal*n1)<0.01,:); 
n2 = mex_ransac_manhattan_2nd(normal_sub',n1); % secndoly get other two axes (should be x and y-axis)
M = [reshape(n2,3,2),n1]';
if M(3,3) < 0
    M = -M;    
end

end