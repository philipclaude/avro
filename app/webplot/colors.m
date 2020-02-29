clc;close all;

n = 256;
%c = hot(n);
%load('viridus.mat')
c = viridis_data;

fid = fopen('colors.h','w');

fprintf(fid,'static int ncolormap = %d;\n',n-1);
fprintf(fid,'static float color_map[%d*3] = \n {',n);
for i=1:n
    fprintf(fid,'%g,%g,%g',c(i,1),c(i,2),c(i,3));
    if i<n
        fprintf(fid,',');
    end
    if mod(i,4)==0
        fprintf(fid,'\n');
    end
end

fprintf(fid,'};\n');

fclose(fid);