

TPN = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\datSkels\';
fileTag = 'vg4';

edgeDat = load([TPN fileTag '_edges_1.dat']);
radDat = load([TPN fileTag '_radii_1.dat']);
vertDat = load([TPN fileTag '_vertices_1.dat']);

subsAll = load([TPN 'cid4sub.mat']);
subs = subsAll.sub;

pred = edgeDat(:,1);
vert = vertDat/100;
edges = edgeDat;
%%

clf

downSamp = 1;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;

smallSub = shrinkSub(subs,downSamp);
%smallSub = smallSub(:,flipDim);

fv = subVolFV(smallSub,[],renderProps);
%fv.vertices = fv.vertices * downSamp;
fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;

[p] = renderFV(fv,[0 0 1],.4);
view([0 0])
axis off
hold on

%scatter3(subs(:,1),subs(:,2),subs(:,3),'.','r')
%scatter3(vert(:,1),vert(:,2),vert(:,3),radDat,'.','g')

plot3([vert(edges(:,1),1) vert(edges(:,2),1)]',[vert(edges(:,1),2) vert(edges(:,2),2)]',...
    [vert(edges(:,1),3) vert(edges(:,2),3)]','w')

hold off






