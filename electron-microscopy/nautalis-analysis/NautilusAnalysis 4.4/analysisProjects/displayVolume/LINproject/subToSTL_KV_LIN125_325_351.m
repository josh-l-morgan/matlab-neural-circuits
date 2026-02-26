clear all

%load('MPN.mat')
MPN = 'D:\LGNs1\Export\export_KV_LIN_2018+09+20a\';

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

cellList = [251 270 ];

col = [1 0 0; 0 1 0; 0 0 1];

renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end


downSamp = 3;

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










