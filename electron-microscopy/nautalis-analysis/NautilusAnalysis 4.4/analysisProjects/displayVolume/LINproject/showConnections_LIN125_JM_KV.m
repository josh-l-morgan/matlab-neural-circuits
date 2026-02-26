
clear all

rotDir = 'C:\Users\jlmorgan\Documents\LIN\images\JM+KV_04\';


%% Load data

MPN = 'D:\LGNs1\Export\export_KV_LIN_morph_2019+3+7K\';
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
obI2 = obI;
dsObj2 = dsObj;


load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir;
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end


    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])


disp('showing')
dSamp = [8 8 4];

colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);
viewProps.dim =3;
viewProps.perspective = 0;


%% color relationship to 124

showCells = [125];
col = [1 0 0];

showCellNames = cat(2,num2cell(showCells));

cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

fsize = double(max(cat(1,dsObj.subs),[],1))+100;
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [0 0 0; fsize];

%% Display Variables

viewProps.maxScaleFactor = .4;
viewProps.sumScaleFactor = .7;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;

viewProps.keepRat = 1;
viewProps.contrast = 2;
viewProps.gamma = .3;
viewProps.dilate = 1;

viewProps2 = viewProps;
viewProps2.cellId = {'125'}
viewProps2.obI = obI2;
viewProps2.dsObj = dsObj2;


degs = [0];
for d = 1:length(degs);
    sprintf('rendering angle %d (%d deg) of %d',d,degs(d),length(degs))
    viewProps.degRot = degs(d);
    
    
    %% Display Cells
    
    %I_topSum = showCellsAndMore(viewProps);
    I_topSum = stereoCellsAndMoreFull_PS(viewProps);
    I_topSum = tweakI(I_topSum,viewProps);
    image(uint8(I_topSum))
    
    I_topSum2 = stereoCellsAndMoreFull_PS(viewProps2);
    I_topSum2 = tweakI(I_topSum2,viewProps);
    image(uint8(I_topSum2))
    
    
    colCom(:,:,1) = sum(I_topSum,3);
    colCom(:,:,2) = sum(I_topSum2,3);
    colCom(:,:,3) = sum(I_topSum2,3);
    
        image(uint8(colCom))

    if 0
        %     rotFold = sprintf('%striad_%04.0f\\',rotDir,useTri);
        iNam = sprintf('tri%04.0f_rot%04.0f.png',0,degs(d));
        if ~exist(rotDir,'dir'),mkdir(rotDir),end
        writeImage = uint8(colCom);%(1000:2500,500:2000,:));
        imwrite(uint8(writeImage),[rotDir iNam])
    end
    
    
end %end rotation
