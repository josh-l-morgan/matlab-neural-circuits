
% clear all



%% Load data
load('MPN.mat')
%MPN = 'D:\KarlsRetina\HxQ\export_2018+03+29c\'
if ~exist('MPN','var')
    MPN = GetMyDir;
end

rotDir = [WPN 'rotDir\rot1'];

synDir = [WPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end

disp('showing')
dSamp =  (obI.em.vRes )./1000/2 .* [1.25 1.25 1.1];
dSamp =  (obI.em.vRes )./1000*8;
dSamp = (obI.em.res )./1000 ./ obI.em.dsRes;
%dSamp = (obI.em.vRes./obI.em.dsRes/1000);

%dSamp = [4 4 40]/[8 8 1] .* [.2 .2 .2] 
%dSamp = obI.em.dsRes./obI.em.res;
%dSamp = 1./dSamp


colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);
viewProps.dim =3;
viewProps.perspective = 0;

%% test
if 1
anchors = double(obI.colStruc.anchors);
sub = cat(1,dsObj(:).subs);
scatter3(sub(:,1),sub(:,2),sub(:,3),'.','b')
hold on
subs2 = anchors;
subs2 = subs2 .* repmat(dSamp,[size(subs2,1) 1]);
scatter3(subs2(:,1),subs2(:,2),subs2(:,3),400,'o','filled','r')
hold off
end


%% Set Variables
targCell = [6];
targCells = targCell;


%% n popsEach possible pairing between RGCs and TCs in the subnetworks of cell A and B was classified according to whether the RGC and TC had the same seed cell association or a different seed cell association.
allEdges = obI.nameProps.edges;
preTarg = preTo(allEdges,targCell);
postTarg = postTo(allEdges,targCell);
%%Reference lists
showCell = [targCell preTarg(:,1)' postTarg(:,1)'];
%showCell = [targCell];

isPost = preTarg(:,1);
isPre = postTarg(:,1);

%% spreadsheet positions
%synPos = getSynPos(1);
%fromRGC_dat = dsAnchors(synPos.postRGCPos,obI,[2 1 3]);


%% motifs

mot = getMotifs(obI,dSamp);
%synStruct = getSynMat(obI);
%anaMot = analyzeMotifs(obI);

%% syns
syns = mot.syn.syns;
fromTarg = syns(:,1) == targCell;
toTarg = syns(:,2) == targCell;

toTargPos = cat(1,mot.syn.synPos{toTarg,:});
fromTargPos = cat(1,mot.syn.synPos{fromTarg,:});


%% triads
if 0
    tri = mot.tri;
    triCell = tri.triCell;
    triCell = triCell(triCell(:,1) == targCell,:);
    triTarg = triCell(:,3);
    triTargRGC = intersect(triTarg,rgcs);
    triTargTCR = intersect(triTarg,tcrs);
    triTargLIN = intersect(triTarg,lins);
    
    
    prim125 = find(tri.triCell(:,1)==125);
    useTri = prim125;
    %for u = 1:length(prim125);%1:size(tri.triCell,1);
    
    triPoints = drawTriads(tri,useTri);
    
    %synCol = cat(1,triPoints.lineCol,triPoints.ballCol);
    % synRad = cat(1,triPoints.lineRad,triPoints.ballRad);
    % synType = cat(1,triPoints.lineRad*0+1,triPoints.ballRad*0+2);
    % useSyn = cat(2,triPoints.lineGroup, triPoints.ballGroup);
    primCell = unique(triPoints.cellGroup(:,1));
    secCell = unique(triPoints.cellGroup(:,2));
    tertCell = unique(triPoints.cellGroup(:,3));
    
end

%% diads
if 0
    diads = mot.di.diCell;
    diads = diads(diads(:,1) == 125,:);
    diTarg = diads(:,3);
end
%% color relationship to 124


group{1} = targCell;
group{2} = [isPost];% setdiff(synTargLIN,125);
group{3} = [isPre];% setdiff(synTargTCR,125);
group{4} = [];% setdiff(unique([diTargTCR; triTargTCR]),125);
groupCol = [1 1 1; 0 1 0; 1 0 0; .5 1 0]*2;

showCells = [];
col = [];
for g = 1:length(groupCol);
    showCells = cat(1,showCells,group{g}(:));
    col = cat(1,col,repmat(groupCol(g,:),[length(group{g}) 1]));
end

showCellNames = cat(2,num2cell(showCells));

cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

fsize = double(max(cat(1,dsObj.subs),[],1))+100;
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [0 0 0; fsize];

%% Color Syns

useSyn{1} = [fromTargPos];%fromRGC_dat;
%useSyn{2} = synPos125fromRGC;
useSyn{2} = [toTargPos];%%fromLIN_dat;
useSyn{3} = [];%synPos125fromUNK;
useSyn{4} = [];%synPos125toTCR;
useSyn{5} = [];%synPos125toLIN;
useSyn{6} = [];%synPos125toUNK;%[];%synPos125fromLIN;


synType = {     'ring'      'bar'       'bar'       'disk'      'ball'        'disk'};
synCol = [      0 1 1;      0 1 1;      0 1 1;      0 1 1;          0 1 1;      0 1 1;];
var1 = [        10          10          64          64              64          64]*2;
var2 = [        5           5           16           16           16           16]


%% Display Variables

viewProps.maxScaleFactor = .2;
viewProps.sumScaleFactor = 1;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;


viewProps.keepRat = .5;
viewProps.contrast = 1.5;
viewProps.gamma = .4;
viewProps.dilate = 0;


degs = [1:1:360];
for d = 1:length(degs);
    sprintf('rendering angle %d (%d deg) of %d',d,degs(d),length(degs))
    viewProps.degRot = degs(d);
    
    
    %% Display Cells
    
    %I_topSum = showCellsAndMore(viewProps);
    I_topSum = stereoCellsAndMoreFull_PS(viewProps);
    I_topSum = tweakI(I_topSum,viewProps);
    %I_topSum = I_topSum*1000;
    
    image(uint8(I_topSum*1))
    
    %% Draw synapses
    
    dim = viewProps.dim;
    if dim == 1
        dims = [3 2];
    elseif dim == 2
        dims = [3 1];
    elseif dim == 3
        dims = [1 2];
    end
    
    %% Draw syn from track sheet
    %synPos = getSynPos(1);
    %     useSyn{1} = synPos.allPos;
    %     useSyn{2} = synPos.prePos;
    %     useSyn{3} = synPos.postRGCPos;
    %     useSyn{4} = synPos.postUnkPos;
    %     useSyn{5} = [11666 13942 3417];
    %     synCol= [ 0 0 0; 2 0 0; 0 1 0; 0 0 3; 0 0 0; 0 0 0]* 1 ;
    %     %synCol= [ 0 0 0; 0 1 1; 1 0 1; 0 0 0; 0 0 0; 0 0 0] ;
    %
    
    showSynTrack  = I_topSum *0;
    
    for p = 1:length(useSyn)
        anchors = double(useSyn{p});
        if sum(anchors(:))
            useAnchors = anchors;
            useAnchors = useAnchors(sum(useAnchors <= 1,2)==0,:);
            ballSum = jmkern(synType{p},var1(p),var2(p));
            viewProps.points = useAnchors;
            useAnchors = stereoPoints_PS(viewProps);
            
            Isize = [viewProps.viewWindow(2,dims(1)) - viewProps.viewWindow(1,dims(1)) ....
                viewProps.viewWindow(2,dims(2)) - viewProps.viewWindow(1,dims(2))];
            
            %         anchorInd = sub2ind(Isize, useAnchors(:,dims(1)) - viewProps.viewWindow(1,dims(1)),...
            %             useAnchors(:,dims(2)) - viewProps.viewWindow(1,dims(1)));
            anchorInd = sub2ind(Isize, useAnchors(:,dims(1)) , useAnchors(:,dims(2)) );
            
            uAnchInd = unique(anchorInd);
            if length(uAnchInd)>1
                valAnch = histc(anchorInd,uAnchInd);
            else
                valAnch = 1;
            end
            
            synImage = zeros(Isize);
            synImage(uAnchInd) = valAnch;
            synImage = convn(synImage,ballSum,'same');
            
            for c = 1:3
                showSynTrack(:,:,c) =  showSynTrack(:,:,c) + synImage * synCol(p,c);
            end
        end
        %         image(uint8(showSynTrack*256));
        %         pause
    end
    
    
    %showSynTrack = outline(showSynTrack,6);
    maskSyn = repmat(sum(showSynTrack,3),[1 1 3]);
    
    showSynTrack = showSynTrack * 256;
    
    showSynTrack(~maskSyn) = I_topSum(~maskSyn);
    % showSynTrack(~maskSyn) = 256 - I_topSum(~maskSyn) * 3;
    %showSynTrack(~maskSyn) = I_topSum(~maskSyn) .* 0 + 255 - I_topSum(~maskSyn)*3;
    %showSynTrack = showSynTrack .* 0 + 255 - showSynTrack;
    
    %showSynTrack = showSynTrack * 256 + I_topSum;
    
    %showSynTrack = 255-showSynTrack;
    
    image(uint8(showSynTrack)),pause(.1)
    
    %     hold on
    %     scatter(useAnchors(:,dims(2)) , useAnchors(:,dims(1)), 'T','w')
    %     pause
    %     hold off
    
    
    if 1
        %     rotFold = sprintf('%striad_%04.0f\\',rotDir,useTri);
        iNam = sprintf('tri%04.0f_rot%04.0f.png',0,degs(d));
        if ~exist(rotDir,'dir'),mkdir(rotDir),end
        writeImage = uint8(showSynTrack);%(1000:2500,500:2000,:));
        imwrite(uint8(writeImage),[rotDir iNam])
    end
    
    
end %end rotation
