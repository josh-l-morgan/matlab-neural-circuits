clear all

load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


obMovDir = [WPN '\movies\IxQ\cid4_VG4_B_VGBip\'];
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end


region =1;
flipDim = [1 3 2 ];

markersOn = 0;
synSegOn = 0;
shouldWrite = 0;
onlyTest = 0;
rotate = 0;
showScale = 0;
fullContext = 0;
randCol = 0;

az = 0:2:360;



downSamp = 1;
renderProps.smooth = 0;
renderProps.resize = 2;
renderProps.smoothPatch = 0;

%% Get cell info
allEdges = obI.nameProps.edges;
allIds = unique(obI.nameProps.cellNum);
allSegs = 1:length(dsObj);
%allCells = obI.cell.name;

sourceCol = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1;1 .4 0; 0 .4 1];
sourceName = obI.fuse.exportDir;
source = obI.fuse.obSource;
alph = source*0+1;
col = sourceCol(source,:);


%% Pick cells
% mot = getMotifs(obI);
% allCells = mot.cel.cells;
% rgcs = mot.cel.types.rgcs;
% tcrs = mot.cel.types.tcrs;
% lins = mot.cel.types.lins;
%
isLocal = [allIds ];%[9079 353 468 221 362 5119];%[ 9078 805 399 426 436 423 9003 434 819 349 492 9193];% 221 353 402]; %% shaft synapses1
group1 = [find(source == 1)];%bpcON;
group2 = [find(source == 2)];
group3 = [find(source == 3)];%[allIds(12)];
group4 = [];%
group5 = [];%setdiff(Post(:,1),[seedList 0])';
group6 = [];%setdiff(allIds,[group1 group2 group3 group4 group5]);


%% Draw cells
% cellList = 125;
aL = lightangle(0,45) ;
trackCells = [];
cellsShown = {};
cellsNotShown = {};
for i = 1:length(dsObj)
    sub = dsObj(i).subs;
    if ~isempty(sub)
        if exist('crop','var')
            if (fullContext == 0) | (i>1)
                useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                    (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                    (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
                sub = sub(useSub,:);
            end
        end
        smallSub = shrinkSub(sub,downSamp);
        smallSub = smallSub(:,flipDim);
        if onlyTest
            smallSub = smallSub(1:400,:);
        end
        %trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            %disp(sprintf('no points on %d',cellList(i)))
            cellsNotShown = cat(1,cellsNotShown,{cellList(i)});
        else
            
            fv = subVolFV(smallSub,[],renderProps);
            [p] = renderFV(fv,col(i,:),alph(i));
            view([0 0])
            axis off
            pause(.01)
            hold on
            %             fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
            %             fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
            %STLWRITE(FILE, FACES, VERTICES)
            %stlwrite(fileNameSTL,fv.faces,fv.vertices);
            %vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
            %toc
            %     cellDat(i).subs = sub;
            %     cellDat(i).fv = fv;
        end
    end
    % disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end
%aL = lightangle(0,45) ;

%trackCells = trackCells(trackCells(:,2)>0,:)


%% movie
tag = 'testMove';
el = zeros(length(az),1);

savefig([obMovDir tag '.fig'])

%
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])


imageName = sprintf('%sspringRun_%s%05.0f_01h.png',obMovDir,tag,i);

% [az el] = view
startAngle =[-88.6 -88.7];
view(startAngle)
lightangle(aL,19.8, -22.8)
pause(.01)
axis off


if 0
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
end


if ~rotate
    view(startAngle)
    lightangle(aL,startAngle(1),startAngle(2)+10)
else
    
    %%
    
    ax = gca;
    camPos = ax.CameraPosition;
    sprintf('%.0f %.0f %.0f',camPos(1),camPos(2),camPos(3))
    camTarg = ax.CameraTarget;
    sprintf('%.0f %.0f %.0f',camTarg(1),camTarg(2),camTarg(3))
    viewDistXY = sqrt((camPos(1)-camTarg(1)).^2 + (camPos(2)-camTarg(2)).^2);
    viewDistXY = 10000./ downSamp;
    camPos = [472 2451 3142]./ downSamp;
    camTarg = [472 2451 3142] ./ downSamp;
    ax.CameraTargetMode = 'manual';
    ax.CameraPositionMode = 'manual';
    ax.CameraPosition = camPos;
    ax.CameraTarget = camTarg;
    campos('manual');
    %ax.Projection = 'orthographic';
    ax.Projection = 'perspective';
    ax.CameraTargetMode
    pause(1)
    for i = 1:length(az);
        i
        if exist('scaleBar'), delete(scaleBar); end
        
        
        viewPosX = sin(deg2rad(az(i)))*viewDistXY;
        viewPosY = cos(deg2rad(az(i)))*viewDistXY;
        ax.CameraPosition = [camTarg(1)+viewPosX camTarg(2) + viewPosY camPos(3)];
        aL.Position = ax.CameraPosition + [ 0 0 5000];
        
        %         view([az(i)+startAngle(1) startAngle(2)])
        %         lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
        %
        %%rotate scalebar
        mz = makehgtform('zrotate',(az(i)+startAngle(1))*2*pi/360);
        %mx = makehgtform('yrotate',0);
        m = mz;
        
        
        pause(.01)
        axis off
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 512 512])
        %runSprings(springDat,allResults{1})
        set(gcf, 'InvertHardCopy', 'off');
        imageName = sprintf('%srot_%05.0f.png',obMovDir,i);
        if shouldWrite
            print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
        end
        
    end
end


