function[] = showRads3D(arbor);







pos = arbor.nodes.pos;
pos = pos - repmat(arbor.skel.offset,[size(pos,1) 1])+1;
rad = arbor.nodes.rad;

subs = arbor.vox.subs;
subs = subs - repmat(arbor.skel.offset,[size(subs,1) 1]) + 1;



renderCon(subs,[],[0 0 1],.2);





[x,y,z] = sphere(31);
fv = surf2patch(x,y,z);
vert = fv.vertices ;
face = fv.faces;

vSize = size(vert,1);
fSize = size(face,1);
pSize = size(pos,1);

vPos = zeros(vSize*pSize,3);
faces = zeros(fSize*pSize,4);

for i = 1:pSize
    newFace = face + (i-1) * vSize;
    newVert = vert * rad(i) + repmat(pos(i,:),[vSize 1]);
    vPos((i-1)*vSize + 1 : (i*vSize),:) = newVert;
    faces((i-1)*fSize + 1 : (i*fSize),:) = newFace;
end

fv.vertices = vPos;
fv.faces = faces;



%%display

vertNum = size(fv.vertices,1);
edgeCol = [.5 .5 .5];
faceCol = [1 0 0];

fv.CDataMapping = 'direct';%repmat([1 1 1],[size(fv.vertices,1) 1]);
fv.FaceVertexCData = faceCol;%repmat(faceCol,[size(fv.vertices,1) 1]);;
fv.FaceColor = 'flat';%[1 1 1];%repmat([1 1 1],[size(fv.vertices,1) 1]);%'flat';
fv.EdgeColor = [0 0 0];
fv.LineWidth = 1;

fv.FaceAlpha = 1;
fv.EdgeAlpha = 0;
fv.AlphaData = 'direct';
fv.AlphaDataMapping = 'direct';

fv.FaceLighting = 'gouraud';
fv.AmbientStrength = .5;
fv.DiffuseStrength = .9;
fv.SpecularExponent = 10;
fv.SpecularStrength = .9;
fv.BackFaceLighting = 'lit';

col = [1 0 0];
alph = .4;

p = patch(fv)

view(30,-15);
axis vis3d;
%colormap copper
set(p,'EdgeColor','none');
p.FaceColor = col;

p.FaceLighting = 'gouraud';
p.AmbientStrength = .2;%.4; shadowColor
p.DiffuseStrength = .8;%.1;
p.SpecularStrength = .6;
p.SpecularExponent = 3;
p.BackFaceLighting = 'lit';
p.Clipping = 'off';
p.FaceAlpha = alph;


daspect([1,1,1])
view(3); axis tight
set(gca,'color',[0 0 0])
set(gcf,'color',[0 0 0])
%set(gca,'off')

l = lightangle(145,45) ;














                
                
                