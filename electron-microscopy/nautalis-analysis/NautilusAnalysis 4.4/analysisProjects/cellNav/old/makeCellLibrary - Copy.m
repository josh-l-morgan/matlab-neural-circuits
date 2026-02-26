clear all

load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
load([WPN 'tis.mat'])


obMovDir = [WPN '\movies\IxQ\cid4_VG4_B_VGBip\'];
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end

flipDim = [1 3 2 ];

downSamp = 1;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;


%% Pick cells

cellList = tis.cells.cids;

col = ones(length(cellList),3);
alph = ones(length(cellList),1);


%%
fvDir = [WPN 'fvLibrary\'];
if ~exist(fvDir,'dir'),mkdir(fvDir),end

save([fvDir 'obI.mat'],'obI')
save([fvDir 'tis.mat'],'tis')

%% Draw cells
trackCells = [];
cellsShown = {};
cellsNotShown = {};
for i = 10:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    obName = cellList(i);
    fvFilename = sprintf('%s%d.mat',fvDir,obName);
    if 1%~exist(fvFilename,'file')
        if ~isempty(sub)
            if iscell(obName); obName = obName{1};end
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
        else
            smallSub = [];
        end
        %trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            fv.vertices = [];
            fv.faces = [];
        else
            fv = subVolFV(smallSub,[],renderProps);
        end
        fvFilename = sprintf('%s%d.mat',fvDir,obName);
        save(fvFilename,'fv');
    end
    
end



