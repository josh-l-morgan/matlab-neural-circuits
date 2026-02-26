function[fv] = meshShape(type,s1,s2,s3);

if ~exist('type','var')
    type = 'cone';
end

if ~exist('s1','var'),    s1 = 1;   end
if ~exist('s2','var'),    s2 = 1;   end
if ~exist('s3','var'),    s3 = 1;   end


if strcmp(lower(type),'square')
    vert =  [-s1 -s1; -s1 s1; s1 s1; s1 -s1];
   %edge = [1 2; 1 4; 3 2; 3 4];    
    face = [ 1 2 3; 1 3 4];
end

if strcmp(lower(type),'sphere')
   [x,y,z] = sphere;
    surf(x,y,z)
end

if strcmp(lower(type),'cone')
    
   a = [0:.1:(2*pi)];
    
   y = sin(a) * s1;
   x = cos(a) * s1;
   c1 = [0 0 s2];
   c2 = [0 0 0];
   circNum = length(y);
   vert = [c1; c2 ; [x' y' zeros(circNum,1)]+s2];
   
   circ1 = 3:size(vert,1);
   circ2  = circshift(circ1,1);
   disk = [circ1' circ2' ones(circNum,1)];
   cone = [circ1' circ2' ones(circNum,1)+1];
   face = cat(1,cone,disk);
   
end


clear fv
vertNum = size(vert,1);
edgeCol = [0 0 0];
faceCol = [1 0 0];

fv.vertices = vert;
fv.faces = face;

fv.CDataMapping = 'direct';%repmat([1 1 1],[size(fv.vertices,1) 1]);
fv.FaceVertexCData = faceCol;%repmat(faceCol,[size(fv.vertices,1) 1]);;
fv.FaceColor = 'flat';%[1 1 1];%repmat([1 1 1],[size(fv.vertices,1) 1]);%'flat';
fv.EdgeColor = [0 0 0];
fv.LineWidth = 1;

fv.FaceAlpha = 1;
fv.EdgeAlpha = 0;
fv.AlphaData = 'direct';
fv.AlphaDataMapping = 'direct';

fv.FaceLighting = 'flat';
fv.AmbientStrength = .5;
fv.DiffuseStrength = .5;
fv.SpecularExponent = 1;
fv.BackFaceLighting = 'unlit';
% 
% patch(fv)
% cl = camlight(45,180);
