


OPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\exports\exportJosh_14+08+12_Mip3_mat\obFiles\'
testOb = imread()
viewWindow = [1000 1000 1000; 5000 5000 5000];

filename = [OPN 'dSamp4_108.obj']

fid = fopen(filename);
obText = fileread(filename);
fclose(fid)


for i=1:size(v,1)
fprintf(fid,'v %f %f %f\n',v(i,1),v(i,2),v(i,3));
end

fprintf(fid,'g %s\n',objectname);

for i=1:size(f,1);
fprintf(fid,'f %d %d %d\n',f(i,1),f(i,2),f(i,3));
end
fprintf(fid,'g\n');

fclose(fid);




p = patch(fv);
view(30,-15);
axis vis3d;
colormap copper
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud