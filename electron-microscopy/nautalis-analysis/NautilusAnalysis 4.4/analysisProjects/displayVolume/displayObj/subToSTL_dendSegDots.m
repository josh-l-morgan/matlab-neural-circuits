clear all

MPN = GetMyDir
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


downSamp = 1;
renderProps.renderOb  = 0;
renderProps.resize = 4; %4
%renderProps.smooth = 1; %1
renderProps.smoothPatch = 0; %3
ballSize = 3;
ballAlf = .65;

%% Get anchors
obNames = obI.colStruc.names;


%obI.em.dsRes = [4 6 30];
anchors = obI.colStruc.anchors;
dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;

anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);
anchors = round(anchors);
%anchors(anchors<1) = 1;



clear isDend isSyn
for i = 1:length(obNames)
    isDend(i) = sum(regexp(obNames{i},'dend'))>0;
    isSyn(i) = sum(regexp(obNames{i},'Segment'))>0;
    isCap(i) = sum(regexp(obNames{i},'cap'))>0;
    isAxon(i) = sum(regexp(obNames{i},'axon'))>0;
    isHighlight(i) = sum(regexp(obNames{i},'highlight'))>0;
end

synAnc = (anchors(isSyn,:));
synAnc = synAnc(sum(synAnc,2)>2,:);
synAnc = synAnc(end:-1:1,:);



%synAnc = getSynAnchors(obI,125, 385);

synBallz = ballz(synAnc,ballSize);
colMap = hsv(256);
numBallz = length(synBallz);

ballzCol = colMap(randperm(256,numBallz),:);

ballzCol = [rand(numBallz,1)*.6  ones(numBallz,1) rand(numBallz,1)*.8];
%ballzCol = repmat([0 1 0],[length(synBallz) 1]);
ballzAlf = repmat(ballAlf,[length(synBallz) 1]);


linked = [];
%cellList = [125 385 ]
cellList = {'dend' 'axon' 'cap' 'highlight' synBallz{:}};
col = [1 1 1; 1 0 0; 0 0 1; 1 1 0; ballzCol];
cellAlpha = [1; 1; 1; .4; ballzAlf];
cellSmooth = [1; 1; 0; 0; ones(length(synBallz),1)*1];

%
% cellList = { cat(1,synBallz{:})};
% col = [ 1 0 0];
% cellAlpha = [ .6]



%%

renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end




target = [120 110 110]*8; % ( X Y Z)
rad = 300;
%crop = [target- rad; target + rad];
%crop     = [ 660   880   680;           1100   1100  1000]; %[ z x y] ; %faciculation
%crop     = [ 760   880   580;           1100   1200  1000]; %[ z x y]

clf
% cellList = 125;
l = lightangle(0,45) ;

for i = 1:length(cellList)
    
    subCell = names2Subs(obI,dsObj,cellList(i));
        sub = subCell{1};

    if ischar(cellList{i})
        if sum(regexp(cellList{i},'highlight'))
            blockSize = [51 51 1];
            sub = highlightSubs(sub,blockSize);
        end
    end
    obName = cellList{i};
    
    if ~ischar(obName)
        if length(obName>1)
            obName = sprintf('object_%d',i);
        end
    end
    
    if iscell(obName); obName = obName{1};end
    if exist('crop','var')
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
        sub = sub(useSub,:);
        
    end
    smallSub = shrinkSub(sub,downSamp);
    %    smallSub = growSub(sub,downSamp);
    smallSub = smallSub(:,[2 3 1]);
    
    
    tic
    if isempty(smallSub)
        disp(sprintf('no points on %d',cellList(i)))
    else
        
        renderProps.smooth = cellSmooth(i);
        fv = subVolFV(smallSub,[],renderProps);
        [p] = renderFV(fv,col(i,:),cellAlpha(i));
        
        view([0 0])
        set(gca,'Zdir','reverse')
        axis off
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
    %disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end

hold off


%% movie
tag = 'dsDend01';
frames = 180;
el = ones(frames,1) * 0;
az = 0:(360/frames):360;%(360-360/frames);
obMovDir = 'Z:\joshm\LGNs1\Analysis\movies\MAC_glom1j\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
savefig([obMovDir tag '.fig'])

%
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])
shouldWrite = 0;
%l = lightangle(0,45) ;



while 1
    for i = 1:frames;
        
        view([az(i) el(i)])
        axis off
        set(gca,'Zdir','reverse')
        lightangle(l,az(i)+20, -30)
        pause(.01)
        
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
        
        %runSprings(springDat,allResults{1})
        set(gcf, 'InvertHardCopy', 'off');
        imageName = sprintf('%sspringRun_%s%05.0f.png',obMovDir,tag,i);
        %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
        
        if shouldWrite
            print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
        end
        
        
    end
    if shouldWrite, break,end
end






