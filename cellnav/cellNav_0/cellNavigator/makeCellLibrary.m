function  makeCellLibrary

global glob tis tisDat
tempFig = figure;
tempAx = gca(tempFig);

MPN = glob.NA.MPN;
WPN = glob.NA.WPN;

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
%load([WPN 'tis.mat'])
tis = makeTis;


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

delete(fvDir)
mkdir(fvDir)

save([fvDir 'obI.mat'],'obI');
save([fvDir 'tis.mat'],'tis');

%% Draw cells
trackCells = [];
cellsShown = {};
cellsNotShown = {};
for i = 1:length(cellList)
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
        
        %% Scale
        fv.vertices = fv.vertices * .1;
        
        fvFilename = sprintf('%s%d.mat',fvDir,obName);
        save(fvFilename,'fv');
        cla
        patch(tempAx,fv)
        pause(.01)
    end
    
end

%% Run SMs

smDir = [WPN 'SMs\'];
libDir = [WPN 'fvLibrary\'];
dSMs = dir([smDir '*.mat']);
nams = { dSMs.name};


for i = 1:length(nams)
    clf
    load([smDir nams{i}])
    %showNepBones(sm.nep);
    %showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.nodeRad,[0 1 0],.4)
    if isfield(sm,'nep')
        if ~isempty(sm.nep)
            fv = showRadSurf_cnv(sm.nep.pos,sm.nep.edges,sm.nep.nodeRad*0+.1,[1 1 1],1,tempFig);
            pause(.01)
            filename = sprintf('%sskelFV_%d.mat',libDir,sm.cid);
            save(filename,'fv')
            
            filename = sprintf('%snep_%d.mat',libDir,sm.cid);
            nep = sm.nep;
            save(filename,'nep')
        end
    end
end
close(tempFig)


renderRefs
makeTisDat



















