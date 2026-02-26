clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

cellList = [unique([9017  527 9018 156 527 397 9013 394  527 394 ...
    9014 9015 156 9015 9016 395 396 9019 9020 9021]) 125];


[postTarg] = getList_125rand;
postAx = [447:456 10034:10036 160 303 184];
cellList = [125 postAx postTarg'];

allEdges = obI.nameProps.edges;
postTarg = preTo(allEdges,125);
preTarg = postTo(allEdges,125);
cellList = [preTarg(:,1)' postTarg(:,1)' 125]
cellList = {125 'frag ax 125'};
cellList = {  125};%{ 125 325 351}


% 
% colMap = hsv(256);
% colNum = length(cellList)-1;
% rainBow = colMap(floor([1:colNum] * 255/(colNum)),:);
% col = [rainBow; [1 1 1]]
col = [1 1 1;repmat([1 0 0],[length(postAx) 1]);repmat([0 1 0],[length(postTarg) 1])]
col = [1 0 0; 0 1 0; 0 0 1];



renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end



downSamp = 3;

target = [908 490 555.5]*2; % ( Y X Z)
rad = 300;
%crop = [target- rad; target + rad];
%crop     = [2000 400 600;          2500 600 700];
clf
for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    obName = cellList(i);
    if iscell(obName); obName = obName{1};end
    if exist('crop','var')
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
        sub = sub(useSub,:);
        
    end
    smallSub = shrinkSub(sub,downSamp);
    tic
    if isempty(smallSub)
        disp(sprintf('no points on %d',cellList{i}))
    else
      fv = subVolFV(smallSub,[],renderOb);
        %fv.vertices =  fv.vertices * obI.em.dsRes(1) *downSamp;
        % renderFV(fv,[1 0 0]);
        
        %l = lightangle(145,45) ;
        fv.FaceColor = col(i,:);
        fv.EdgeColor = 'none';
        fv.FaceLighting = 'gouraud';
        fv.AmbientStrength = .5;%.4;
        fv.DiffuseStrength = .9;%.1;
        fv.SpecularStrength = 1;
        fv.SpecularExponent = 5;
        fv.BackFaceLighting = 'lit';
        fv.FaceAlpha = 1;
        p = patch(fv)
        
         view(30,-15);
        axis vis3d;
        daspect([1,1,1])
        view(3); axis tight
        camlight
        lighting gouraud
        
        set(gca,'color',[1 1 1])
        set(gcf,'color',[1 1 1])
                
        pause(.01)
        tag = 'cell3'
    %     cellDat(i).subs = sub;
    %     cellDat(i).fv = fv;
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList{i},i,length(cellList)));
    pause(.01)
end
camlight
hold off










