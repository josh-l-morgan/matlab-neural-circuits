
%    MPN = GetMyDir


load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])
cellList = obI.cell.name(obI.cell.isCell>0);
%cellList =

%% Skeletonize with shortest paths
cellList = obI.cell.name(obI.cell.isCell>0);

disp(sprintf('Cell List = %s',num2str(cellList)));
TPN = [MPN 'skel\'];
TPNview = [TPN 'view\'];
TPNmat = [TPN 'mat\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
if ~exist(TPNview,'dir'),mkdir(TPNview);end
if ~exist(TPNmat,'dir'),mkdir(TPNmat);end
morphDir = [MPN 'morph\'];
if ~exist(morphDir,'dir'),mkdir(morphDir),end

morphConDir = [MPN 'morphCon\']
if ~exist(morphConDir,'dir'),mkdir(morphConDir),end


%% Get Preferences
seedList = [108 201 170];
preNames = [1028 1032 2033 2034 2035 2003 2007]


useList = obI2cellList_seedInput(obI,seedList);
conTo = makeConTo(obI,seedList);
seedPref = seedPreferences(seedList,useList);

cellList = useList.postList;
cellCon =  seedPref.cellPrefNoExclusion;
cellPref = seedPref.cellPrefNoExclusion(1,:)./sum(seedPref.cellPrefNoExclusion(1:2,:));
isPref = sum(seedPref.cellPrefNoExclusion(1:2,:));

%%
c = 0;
clear recRes latRat

for i = 1:length(cellList)
    clf
    subplot(2,1,1)

    cellTarg = cellList(i);
    disp(sprintf('skeletonizing cell %d (%d of %d)',cellTarg,i,length(cellList)))
    
    skelFile = sprintf('%s%d.mat',TPNmat,cellTarg);
    viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    
    
    rawObjectSubs = getCellSubs(obI,dsObj,cellTarg);
    rawSeed   = ceil(getSeed(obI,cellTarg));
    
    if (size(rawObjectSubs,1)>1000) & (~isempty(rawSeed))
        pass = 1;
        cellFailed = [];
        
        objectSubs = double(downSampSub(rawObjectSubs,[4 4 4]));
        probableDownSamp = rawSeed./(median(rawObjectSubs)./[4 4 4]);
        seedSub = round(rawSeed./[32 32 16]);
        %cellStruct = subs2arbor(objectSubs,seedSub);
        
        try load(skelFile); %Get cellStruct
        catch err
            pass = 0;
            cellTarg
            err
        end
        
        if pass
            subplot(2,1,1)
            'cell passed'
            cellView = cellStruct.sideViews{1};
            skel = cellStruct.skel;
            image(cellView),pause(.01)
            
            cellArbors.cellName(i) = cellTarg;
            cellArbors.arbors(i) = rmfield(cellStruct.arbor,'vox');
            stopTime = clock;
            
            cellOri = princompSkeletons(cellArbors.arbors(i),seedSub);
%             subplot(2,1,2)
%             indPref = cellCon(:,i);
%             prefColor = [indPref/max(indPref)];
%             bar(indPref,'FaceColor', prefColor)
%             ylim([0 30])
%                      
                        pause(.1)

                fileName = sprintf('%smorphOri_%03.0f.png',morphConDir,cellTarg)
               % saveas(gcf,fileName,'png')
            
            dirCount(:,:,i) = cellOri.dirCount;
            lengthUsed(i) = cellOri.lengthUsed;
            oriRat(:,i) = cellOri.oriRat;
            foundSkel(i) = 1;
            
            c = c+1;
            
            recRes(c,1) = cellTarg;
            recRes(c,2) = cellOri.oriRat(1);
            recRes(c,3) = cellOri.oriRat(2);
             latRes(c,1) = cellTarg;
            latRes(c,2) = cellOri.latRat(1);
            latRes(c,3) = cellOri.latRat(2);
            
        else
            'cell Failed'
            
            dirCount(:,:,i) = [0 0; 0 0; 0 0];
            oriRat(:,i) = [0 0];
            lengthUsed(i) = 0;
            foundSkel(i) = 0;
            
        end
        
        
    end
    
end



%%
subplot(1,1,1)
useCell = (foundSkel>0)  & (isPref>2) & (lengthUsed>300);
scatter(cellPref(useCell), oriRat(1,useCell))

xlim([-.2 1.2])
ylim([0 1])

clf
subplot(3,1,1)
hist(oriRat(1,(useCell & (cellPref==0))))

subplot(3,1,2)
hist(oriRat(1,(useCell & (cellPref>0) & (cellPref<1))))

subplot(3,1,3)
hist(oriRat(1,(useCell & (cellPref==1))))










