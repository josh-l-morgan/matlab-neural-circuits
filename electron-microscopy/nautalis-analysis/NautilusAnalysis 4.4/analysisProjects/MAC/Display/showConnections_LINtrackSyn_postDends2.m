
clear all

rotDir = 'C:\Users\jlmorgan\Documents\MAC\images\postsSlice\';


%% Load data
load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir;
end
    MPN = GetMyDir;

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
preTarg = preTo(allEdges,1007);
postTarg = postTo(allEdges,1007);
%%Reference lists
showCell = [1007 11007];

% isPost = preTarg(:,1);
% isPre = postTarg(:,1);

allPost = unique(allEdges(:,1));
allPre = unique(allEdges(:,2));

otherPost = setdiff(allPost,[1007 11007 0]);


%% spreadsheet positions

synPos = getSynPos(1);

fromRGC_dat = dsAnchors(synPos.postRGCPos,obI,[2 1 3]);
fromLIN_dat = dsAnchors(synPos.postLinPos,obI,[2 1 3]);




%% color relationship to 124


group{1} = [allPost];
group{2} = [];% setdiff(synTargLIN,125);
group{3} = [];% setdiff(synTargTCR,125);
group{4} = [];% setdiff(unique([diTargTCR; triTargTCR]),125);
groupCol = [10 10 10; 0 4 0; 0 0 4; 41 0 0];

showCells = [];
col = [];
for g = 1:size(groupCol,1);
    showCells = cat(1,showCells,num2cell(group{g}(:)));
    col = cat(1,col,repmat(groupCol(g,:),[length(group{g}) 1]));
end

showCellNames = cat(2,num2cell(showCells));
showCellNames = showCells;

cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

fsize = double(max(cat(1,dsObj.subs),[],1))+100;
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [0 0 0; fsize];
%viewProps.viewWindow = [500 1100 200; 1100 1700 800];

viewProps.viewWindow = [500 1000 605; 1300 1700 609];
viewProps.perspective = 0;

4053/9802 * fsize(3)

%% Color Syns
useSyn{1} = [];%fromRGC_dat;
%useSyn{2} = synPos125fromRGC;
useSyn{2} = [];
useSyn{3} = [];%synPos125fromUNK;
useSyn{4} = [];
useSyn{5} = [];%synPos125toLIN;
useSyn{6} = [];


synType = {  'bar'      'bar'       'bar'       'ring'      'ring'      'x'};
synCol = [   0 1 0;    .2 .2 1;    1 0 1;      0 1 0;     .2 .2 1;     1 0 0];
var1 = [    16          24          16          16          24          30];
var2 = [    4           6           4           3           3           30]


%% Display Variables

viewProps.maxScaleFactor = .3;
viewProps.sumScaleFactor = .5;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;

viewProps.keepRat = .3;
viewProps.contrast = .2;
viewProps.gamma = .3;



degs = [0];
for d = 1:length(degs);
    sprintf('rendering angle %d (%d deg) of %d',d,degs(d),length(degs))
    viewProps.degRot = degs(d);
    
    
    %% Display Cells
    
    %I_topSum = showCellsAndMore(viewProps);
    I_topSum = stereoCellsAndMoreFull_PS(viewProps);
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
    
    maskSyn = repmat(sum(showSynTrack,3),[1 1 3]);
    
    showSynTrack = showSynTrack * 256;
    showSynTrack(~maskSyn) = I_topSum(~maskSyn);
   % showSynTrack(~maskSyn) = 256 - I_topSum(~maskSyn) * 3;
    %showSynTrack(~maskSyn) = I_topSum(~maskSyn) .* 0 + 255 - I_topSum(~maskSyn)*3;
    %showSynTrack = showSynTrack .* 0 + 255 - showSynTrack;
    
    %showSynTrack = showSynTrack * 256 + I_topSum;
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
