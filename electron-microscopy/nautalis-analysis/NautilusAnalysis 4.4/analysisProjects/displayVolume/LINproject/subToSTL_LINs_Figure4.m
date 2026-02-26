clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])



downSamp = 2;

region = 1;
markersOn = 0;
shouldWrite = 1;
renderProps.smooth = 1;
%renderProps.resize = 1;
renderProps.smoothPatch = 1;
%renderProps.dilate = 0;
onlyTest = 0;

cellList = {  125 356 208 209 904};


col = [0 0 0; 1 0 1; 0 0 1; 1 0 0; 0 .8 0];



renderOb = 0;
tag = 'testCrop';
obMovDir = 'C:\Users\jlmorgan\Documents\LIN\images\outputLINsSTL_09\';
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end





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
    smallSub = smallSub(:,[ 3 2 1]);
    tic
    if isempty(smallSub)
        disp(sprintf('no points on %d',cellList{i}))
    else
      fv = subVolFV(smallSub,[],renderProps);
        %fv.vertices =  fv.vertices * obI.em.dsRes(1) *downSamp;
        % renderFV(fv,[1 0 0]);
        
        %l = lightangle(145,45) ;
        fv.FaceColor = col(i,:);
        fv.EdgeColor = 'none';
        fv.FaceLighting = 'gouraud';
        fv.AmbientStrength = .5;%.4;
        fv.DiffuseStrength = .9;%.1;
        fv.SpecularStrength = 0;
        fv.SpecularExponent = 5;
        fv.BackFaceLighting = 'lit';
        fv.FaceAlpha = 1;
        p = patch(fv)
        
         view(0,0);
        axis vis3d;
        daspect([1,1,1])
        view(3); axis tight
        aL = camlight
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


view(0,180);
startAngle = [0 180];

 az = 0
 
    for i = 1:length(az);
        i
        view([az(i)+startAngle(1) startAngle(2)])
        lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
        %lightangle(aL,82.4,110)
        
        axis off
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 512 512])
        %runSprings(springDat,allResults{1})
        set(gcf, 'InvertHardCopy', 'off');
        imageName = sprintf('%srot_%05.0f.png',obMovDir,i);
        %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
        %print(gcf, imageName, '-dpng','-opengl','-r72')
        if shouldWrite
            print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
        end
        pause(.01)
    end


