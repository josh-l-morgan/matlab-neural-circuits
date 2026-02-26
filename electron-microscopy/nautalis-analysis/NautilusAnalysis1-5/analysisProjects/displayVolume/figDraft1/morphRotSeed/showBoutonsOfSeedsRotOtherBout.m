

%% Load data
load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir
end

synDir = ['D:\LGNs1\Analysis\seedBouton\morphRotation\test1\'];
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

%% Set Variables
seedCells = [108];
%seedCells = [201];


for tc = 1:length(seedCells)
    clear viewProps
    
    targCells = [seedCells(tc)];
    useList = obI2cellList_seedInput_RGC_TCR(obI,targCells);
    
    
    
    %% Find connectivity
    conTo = makeConTo(obI,targCells);
    
    %% pick cells
    
    %%Reference lists
    
    axnet1 =  [1000:1500];
    axnet2 = [2000:2031 2038:2300];
    axnet3 = [2032:2037];
    oneAx = 2040;
    
    showCell = [targCells ];
    
    %showCell = [targCells conTo(1).rgcList];
    
    showCellNames = cat(2,num2cell(showCell));
    
    viewProps.dim =1;
    col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
    col = col(randperm(size(col,1)),:);
    col(1,:) = 1;
    
    
    %% Pick synapses
    prePop{1} = [conTo(1).rgcList];
    postPop{1} =  [setdiff(useList.postList(:)',targCells)];
    
    %
    % prePop{3} = [];
    % postPop{3} =  [];
    %
    % prePop{1} = [];
    % postPop{1} =  [];
    
    
    cellPicDir = [MPN '\cellPic\'];
    if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end
    
    cellList = obI.cell.name;
    
    
    fsize = double(max(cat(1,dsObj.subs),[],1));
    minVal = double(min(cat(1,dsObj.subs),[],1));
    viewProps.viewWindow = [minVal; fsize];
        viewProps.viewWindow = [1 1 1; fsize];

%         viewProps.viewWindow = [650 1 1300; 2750 fsize(2) fsize(3)];
%         viewProps.viewWindow = [1200 1 1; fsize(1) fsize(2) fsize(3)];
    viewProps.viewWindow = [1200 1 1; 3000 fsize(2) fsize(3)];

    %}
    
    
    %% Display Variables
    
    viewProps.maxScaleFactor = 3;
    viewProps.sumScaleFactor = .1;
    viewProps.obI = obI;
    viewProps.dsObj = dsObj;
    viewProps.col = col;
    viewProps.fsize = fsize;
    viewProps.cellId = showCellNames;
    
    
    
    %% Display Cells
    
%     I_topSum = showCellsAndMore(viewProps);
%     image(uint8(I_topSum*1))
%     
    
    %% Draw boutons
    
    otherList = useList.preList(:)';
    %viewProps.cellId = num2cell([targCells useList.preList]);
    viewProps.cellId = num2cell([targCells otherList]);

    %viewProps.cellId = num2cell([targCells prePop{1}]);
    %viewProps.cellId = num2cell([targCells ]);

    
    col2 = colMap(ceil((1:length(otherList ))*256/length(otherList)),:);
    col2 = cat(1,[1 1 1], col2(randperm(size(col2,1)),:));
    load('.\data\colormaps\col2_ax5.mat')
    % save('.\data\colormaps\col2_ax5.mat','col2')
    
    %col2 = [1 1 1];
    viewProps.col = col2;
    
    
    
    if 0
    I_axons = showCellsAndMore(viewProps);
    image(uint8(I_axons))
    end
    %% Ball
    ballRad = 16;
    ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
    [y x z] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
    ball(:) = dists;
    ball(ball>ballRad) = 0;
    ballSum = sum(ball>0,3);
    ballSum = ballSum * 200/max(ballSum(:));
    [y x z ] =ind2sub(size(ball),find(ball>0));
    ballSub = [y x z] - ballRad - 1;
    
    dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
    
    anchors = double(obI.colStruc.anchors);
    anchors(:,1) = anchors(:,1)*dSamp(1);
    anchors(:,2) = anchors(:,2)*dSamp(2);
    anchors(:,3) = anchors(:,3)*dSamp(3);
    anchors = round(anchors);
    anchors(anchors<1) = 1;
    synapses = obI.nameProps.edges(:,1:3);
    
    
    dim = viewProps.dim;
    if dim == 1
        dims = [3 2];
    elseif dim == 2
        dims = [3 1];
    elseif dim == 3
        dims = [1 2];
    end
    
%    showSyn  = I_topSum *0;
    
    usePre = [];
    for i = 1:length(prePop{1})
        usePre = cat(1,usePre,find(synapses(:,2)== prePop{1}(i)));
    end
    
    usePost = [];
    for i = 1:length(postPop{1})
        usePost = cat(1,usePost,find(synapses(:,1)== postPop{1}(i)));
    end
    
    foundSyn = intersect(usePre,usePost);
    
    useSynID = synapses(foundSyn,3);
    useAnchors = anchors(useSynID,:);
    
    dilateAnchors = zeros(size(useAnchors,1)*size(ballSub,1),3);
    S = size(ballSub,1);
    for a = 1:size(useAnchors,1)
        L = (a-1) * S;
        dilateAnchors(L+1:L+S,1) = ballSub(:,1) + useAnchors(a,1);
        dilateAnchors(L+1:L+S,2) = ballSub(:,2) + useAnchors(a,2);
        dilateAnchors(L+1:L+S,3) = ballSub(:,3) + useAnchors(a,3);
    end
    
    %%Get unique
    dilateAnchors(dilateAnchors<1) = 1;
    maxD = max(dilateAnchors,[],1);
    DAind = sub2ind(maxD,dilateAnchors(:,1),dilateAnchors(:,2),dilateAnchors(:,3));
    uDAind = unique(DAind);
    [y x z] = ind2sub(maxD,uDAind);
    dilateAnchors = [y x z];
    
    viewProps.mask = dilateAnchors;
    viewProps.useMask = ones(length( viewProps.cellId),1);
    viewProps.useMask(1) = 0;
    
    
    if 0
        I_boutons = showCellsAndMore(viewProps);
        image(uint8(I_boutons*3))
        
        %%
        % I = Ist;
        I = I_boutons;
        viewProps.keepRat = 1;
        viewProps.contrast = 1;
        viewProps.gamma = 1;
        viewProps.dilate = 0;
        
        I = tweakI(I,viewProps);
        image(uint8(I))
        
        %%
        %   imwrite(uint8(I),[cellPicDir 'ABmix.png'])
        
    end
    
    
    %% 3d display
       
    if 0
        Ist = showStereo(viewProps);
        image(uint8(Ist) )
        %
    end
    
    if 1
        %%
        viewProps.keepRat = 1;
        viewProps.contrast = 1;
        viewProps.gamma = 1;
        viewProps.dilate = 0;
        viewProps.perspective = 0;

        
        rotDir = 'D:\LGNs1\Analysis\movies\morphBout\rot180d\otherBout\'
        if ~exist(rotDir,'dir'),mkdir(rotDir), end
        
        degs = [360:-1:0];
        %degs = [-1.5 1.5];
        %degs = [0:1:359];
        %degs = [0:3:15];
        degs = 0;
        
        I_stack = showRotationPerspective(viewProps,degs,rotDir);
        %playRot(I_stack);
        %%
    end
    
    
    
    
   return 
    
    
end
    
    