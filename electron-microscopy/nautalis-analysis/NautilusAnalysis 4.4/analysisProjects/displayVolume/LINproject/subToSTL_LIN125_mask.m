clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

region = 1;
downSamp = 1;

allEdges = obI.nameProps.edges;

mot = getMotifs(obI);
rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;
unks = mot.cel.types.unks;

seedList = 125;
Pre = postTo(allEdges,seedList);
Post = preTo(allEdges,seedList);

TCR = setdiff(intersect([Post(:,1)],tcrs),seedList)';
RGC = setdiff(intersect([Pre(:,1)],rgcs),seedList)';
LINout = setdiff(intersect([Post(:,1)],lins),seedList)';
LINin = setdiff(intersect([Pre(:,1)],lins),seedList)';
UNKout = setdiff(intersect([Post(:,1)],unks),seedList)';
UNKin = setdiff(intersect([Pre(:,1)],unks),seedList)';

isLocal = []; 

%% shrink mask
shrink = 4;
lookDist = 5;
ballRad = round(lookDist/shrink/ .2);
lookDist = ballRad * shrink * .2;

subCell = names2Subs(obI,dsObj,125);
maskSub = subCell{1};
smallSub = shrinkSub(maskSub,shrink);

volSize = max(cat(1,dsObj.subs),[],1);
smallVolSize = shrinkSub(volSize,shrink);
smallVol = zeros(smallVolSize);

inds = sub2ind(smallVolSize,smallSub(:,1),smallSub(:,2),smallSub(:,3));
smallVol(inds) = 1;


bSub = ballz(smallSub, ballRad);
bSubs = round(cat(1,bSub{:}));
bSubs(bSubs < 1) = 1;
wall = repmat(smallVolSize,[size(bSubs,1) 1]);
bSubs(bSubs > wall) = wall(bSubs >wall);
bInds = sub2ind(smallVolSize,bSubs(:,1),bSubs(:,2),bSubs(:,3));

smallVol(bInds) = 1;



%%
    


maskDist = 5/.2;

group1 = RGC;
group2 = TCR;
group3 = [LINin LINout];
group4 = UNKin;%[intersect([UNKin];%[227]
% 
% group1 = 156
% group2 = 527
% group3 = 9019;
%
% group1 = group1(1);
% group2 = group2(1);

cellList = [125 group1 group2 group3 group4];
% colNum = length(cellList)-1;
% colMap = hsv(256);
% rainBow = colMap(ceil([1:colNum] * 255/(colNum)),:);
% rainBow = rainBow(randperm(size(rainBow,1)),:);
%col = [rainBow; [1 1 1]]
groupCol = [1 1 1; 0 1 0; .3 .3 1; 1 0 0; 1 1 0];


groupAlph = [1 1 1 1 1];
    
col = [groupCol(1,:);repmat(groupCol(2,:),[length(group1) 1]); ...
    repmat(groupCol(3,:),[length(group2) 1]); repmat(groupCol(4,:),[length(group3) 1]);...
    repmat(groupCol(5,:),[length(group4) 1])];
alph = [groupAlph(1);repmat(groupAlph(2),[length(group1) 1]);repmat(groupAlph(3),[length(group2) 1]);...
    repmat(groupAlph(4),[length(group3) 1]); repmat(groupAlph(5),[length(group4) 1])];
% 
% %col = col + rand(size(col))*.5;
% col(2:5,:) = col(2:5,:) + [0 0 0; .7 0 0; 0 0 .7; .7 0 .7];
% col(1,:) = [1 0 0]
% col(col>1) = 1;

%col = [1 1 1];
%%

renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end





% target = [908 490 555.5]*2; % ( X Y Z)
% rad = 200;
% crop = [target- rad; target + rad];
% crop     = [1600  850 850;          1950 1100 1430]; %[ z x y
% crop     = [1750  900 950;          1920 1100 1250]; %[ z x y
% crop     = [1750  900 950;          1920 1100 1250]; %[ z x y



clf
% cellList = 125;
l = lightangle(0,45) ;
trackCells = [];
for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    
    %% check mask
    sub2 = round(sub/shrink);
    sub2(sub2<1) = 1;
    inds = sub2ind(smallVolSize,sub2(:,1),sub2(:,2),sub2(:,3));
    useSub  = smallVol(inds);
    sub = sub(useSub>0,:);
    
    obName = cellList(i);
    if ~isempty(sub)
        if iscell(obName); obName = obName{1};end
        if exist('crop','var')
            useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
            sub = sub(useSub,:);
            
        end
        smallSub = shrinkSub(sub,downSamp);
        trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            disp(sprintf('no points on %d',cellList(i)))
        else
            
            fv = subVolFV(smallSub,[],0);
            
%             if exist('maskSub','var') & (i>1)
%                 vert = fv.vertices;minDist = zeros(size(sub,1),1);
%                 for s = 1:size(vert,1)
%                     subDist = sqrt((maskSub(:,1) - vert(s,1)).^2 + (maskSub(:,2) -vert(s,2) ).^2 ...
%                         + (maskSub(:,3)- vert(s,3) ).^2);
%                     minDist = min(subDist);
%                 end
%                 fv.FaceVertexAlphaData(s,1) = minDist<=maskDist;
%             else
%                       fv.FaceVertexAlphaData = repmat(alph(i),[size(fv.vertices,1),1 1])*alph(i);
%       
%             end
            
            
            fv.facecolor = col(i,:);
            fv.faceAlpha = alph(i);%'interp';%'interp';
            fv.edgeColor = 'none';
          % patch(fv)
%             alphAll = repmat(alph(i),[size(fv.faces,1),1 1]);
            [p] = renderFV(fv);
            view([0 0])
            axis off
            
            
            set(gca,'CameraViewAngle',7,'projection','perspective')
       
        daspect([1 1 1]);
        axis off
        grid off
        set(gcf,'color',[0 0 0])
               
            
            pause(.01)
            hold on
            fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
            fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
            %STLWRITE(FILE, FACES, VERTICES)
            %stlwrite(fileNameSTL,fv.faces,fv.vertices);
            vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
            toc
            %     cellDat(i).subs = sub;
            %     cellDat(i).fv = fv;
        end
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
    end

hold off
trackCells = trackCells(trackCells(:,2)>0,:)

%% movie
tag = 'testMove';
frames = 360;
el = ones(frames,1) * 30;
az = 1:360/frames:360;
obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_Clear_4\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
savefig([obMovDir tag '.fig'])

%
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])
shouldWrite = 0;


imageName = sprintf('%sspringRun_%s%05.0f_01h.png',obMovDir,tag,i);

if region == 1
    % [az el] = view
    
    view(19.8, -62.8)
    lightangle(l,19.8, -22.8)
    pause(.01)
    axis off
    
else
    
    view([82.4 32.4])
    lightangle(l,82.4,110)
    pause(.01)
    axis off
    
end
return
if shouldWrite
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
end

        
% 
% 
% while 1
%     for i = 94%1:frames;
%         
%         view([az(i) el(i)])
%         lightangle(l,az(i)+10, 50)
%         pause(.01)
%         axis off
%         set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
%         
%         %runSprings(springDat,allResults{1})
%         set(gcf, 'InvertHardCopy', 'off');
%         imageName = sprintf('%sspringRun_%s%05.0f.png',obMovDir,tag,i);
%         %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
%         
%         if shouldWrite
%             print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
%         end
%         
%         
%     end
%     if shouldWrite, break,end
%     break
% end
% 
% 




