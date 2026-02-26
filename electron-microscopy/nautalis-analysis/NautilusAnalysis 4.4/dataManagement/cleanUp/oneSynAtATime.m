
clear all
%% Load data
load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir;
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end

disp('showing')
dSamp = [8 8 4];

colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);
viewProps.dim =3;

%% Set Variables
targCells = [108];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n popsEach possible pairing between RGCs and TCs in the subnetworks of cell A and B was classified according to whether the RGC and TC had the same seed cell association or a different seed cell association.
allEdges = obI.nameProps.edges;
preTarg = preTo(allEdges,125);
postTarg = postTo(allEdges,125);
%%Reference lists
showCell = [125 preTarg(:,1)' postTarg(:,1)'];
showCell = [125];

isPost = preTarg(:,1);
isPre = postTarg(:,1);

%% spreadsheet positions

 synPos = getSynPos(1);
   
fromRGC_dat = dsAnchors(synPos.postRGCPos,obI,[2 1 3]);
fromLIN_dat = dsAnchors(synPos.postLinPos,obI,[2 1 3]);


%% motifs

mot = getMotifs(obI);
synStruct = getSynMat(obI);
%anaMot = analyzeMotifs(obI);

rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;

%% syns
syns = [synStruct.pre synStruct.post];
from125 = synStruct.pre == 125;
to125 = synStruct.post == 125;
toTCR  = from125 &(synStruct.postClass== 2);
toLIN = from125 & (synStruct.postClass == 3);
fromRGC = to125 & (synStruct.preClass == 1);
fromLIN = to125 & (synStruct.preClass == 3);
fromUNK = to125 & (synStruct.preClass == 4);


syns = syns(to125,:);
synTarg = syns(:,1);
synTargRGC = intersect(synTarg,rgcs);
synTargTCR = intersect(synTarg,tcrs);
synTargLIN = intersect(synTarg,lins);

synPos125toTCR = synStruct.synPosDS(toTCR,:);
synPos125toLIN = synStruct.synPosDS(toLIN,:);

synPos125fromRGC = synStruct.synPosDS(fromRGC,:);
synPos125fromLIN = synStruct.synPosDS(fromLIN,:);
synPos125fromUNK = synStruct.synPosDS(fromUNK,:);

%% triads
tri = mot.tri;
triCell = tri.triCell;
triCell = triCell(triCell(:,1) == 125,:);
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



%% diads

diads = mot.di.diCell;
diads = diads(diads(:,1) == 125,:);
diTarg = diads(:,3);
diTarg = setdiff(diTarg,triTarg);

diTargRGC = intersect(diTarg,mot.cel.types.rgcs);
diTargTCR = intersect(diTarg,mot.cel.types.tcrs);
diTargLIN = intersect(diTarg,lins);


%% color relationship to 124

postTCR = intersect(tcrs,isPost);
triTarg = intersect(tertCell,tcrs);
diTarg = setdiff(diTarg,triTarg);


group{1} = 125;
group{2} = [];% setdiff(synTargLIN,125);
group{3} = [];% setdiff(synTargTCR,125);
group{4} = [];% setdiff(unique([diTargTCR; triTargTCR]),125);
groupCol = [3 3 3; 0 5 2; 0 0 1; 1 0 0];

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

useSyn{1} = fromRGC_dat;
%useSyn{2} = synPos125fromRGC;
useSyn{2} = fromLIN_dat;
useSyn{3} = synPos125fromUNK;
useSyn{4} = synPos125toTCR;
useSyn{5} = synPos125toLIN;


synType = [ 1 1 1 2 2 2];
synCol = [2 .5 0 ; 1 0 0; 0 0 0;  0.2 .2 1; 1 0 0];
synRad = [ 10 5 5 7 7];


%% Display Variables

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = .5;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;

viewProps.keepRat = .3;
viewProps.contrast = .1;
viewProps.gamma = .3;


synList = synPos125toTCR;
degs = [0];
p = 1;

    viewProps.degRot = 0;

    I_topSum = stereoCellsAndMoreFull_PS(viewProps);

for d = 1:length(synList);
    sprintf('rendering angle %d (%d deg) of %d',d,p,length(synList))
    
    
    %% Display Cells
    
    %I_topSum = showCellsAndMore(viewProps);
    %I_topSum = tweakI(I_topSum,viewProps);

    %image(uint8(I_topSum*1))
    
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
    
        anchors = synList(451,:);
        if sum(anchors(:))
        useAnchors = anchors;
        
        if synType(p) == 1
            ballSum = ballPic(synRad(p));
        else
            ballSum = ringPic(synRad(p));
        end
        
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

    maskSyn = repmat(sum(showSynTrack,3),[1 1 3]);
    
    showSynTrack = showSynTrack * 256;
    showSynTrack(~maskSyn) = I_topSum(~maskSyn);
    %showSynTrack(~maskSyn) = I_topSum(~maskSyn) .* 0 + 255 - I_topSum(~maskSyn)*3;
    %showSynTrack = showSynTrack .* 0 + 255 - showSynTrack;
    
    %showSynTrack = showSynTrack * 256 + I_topSum;
    image(uint8(showSynTrack)),pause(.1)
    pause(.1)
    if 0
        rotDir = 'C:\Users\jlmorgan\Documents\LIN\images\inputTypesSynPos_05\';
        %     rotFold = sprintf('%striad_%04.0f\\',rotDir,useTri);
        iNam = sprintf('tri%04.0f_rot%04.0f.png',0,degs(d));
        if ~exist(rotDir,'dir'),mkdir(rotDir),end
        writeImage = uint8(showSynTrack);%(1000:2500,500:2000,:));
        imwrite(uint8(writeImage),[rotDir iNam])
    end
    
    
end %end rotation
