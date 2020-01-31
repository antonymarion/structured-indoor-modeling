function simple_write_mesh_ply(A, tri, Path)

[fid,Msg] = fopen(Path,'wt');

%%% write PLY header %%%
Format = 'ascii';
fprintf(fid,'ply\nformat %s 1.0\n',Format);
fprintf(fid,'element vertex %d\n',size(A,1));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'element face %d\n', size(tri,1));
fprintf(fid,'property list uint8 int32 vertex_indices\n');
fprintf(fid,'end_header\n');
B = [3*ones(size(tri,1),1), tri];
% write element dat

fprintf(fid, '%g %g %g\n', A');
fprintf(fid, '%d %d %d %d\n',B');
fclose(fid);
end

