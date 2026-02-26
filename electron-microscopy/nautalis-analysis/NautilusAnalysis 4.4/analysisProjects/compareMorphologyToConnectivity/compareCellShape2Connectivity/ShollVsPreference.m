
%    MPN = GetMyDir


load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])
cellList = obI.cell.name(obI.cell.isCell>0);
%cellList =

%figSize = get(gcf,'Position')
figSize = [20 20 900 1000];
set(gcf,'Position',figSize);


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

shollConDirImage = [MPN 'shollCon\image\']
if ~exist(shollConDirImage,'dir'),mkdir(shollConDirImage),end
shollConDirData = [MPN 'shollCon\data\']
if ~exist(shollConDirData,'dir'),mkdir(shollConDirData),end

%% Get Preferences
seedList = [108 201 170 110];
preNames = [1028 1032 2033 2034 2035 2003 2007]


useList = obI2cellList_all(obI,seedList);
conTo = makeConTo(obI,seedList);
seedPref = seedPreferences(seedList,useList);

cellList = useList.postList;
cellCon =  seedPref.cellPrefNoExclusion;
cellPref = seedPref.cellPrefNoExclusion(1,:)./sum(seedPref.cellPrefNoExclusion(1:2,:));
isPref = sum(seedPref.cellPrefNoExclusion(1:2,:));

%%
clear allDOi
overWrite = 1;
for i = 1:length(cellList)
    clf
    cellTarg = cellList(i);
    fileName = sprintf('%smorphSholl_%03.0f.png',shollConDir,cellTarg)
    
    if ~exist(fileName,'file') | overWrite
        
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
                
                'cell passed'
                cellView = cellStruct.sideViews{1};
                skel = cellStruct.skel;
                image(cellView),pause(.01)
                
                cellArbors.cellName(i) = cellTarg;
                cellArbors.arbors(i) = rmfield(cellStruct.arbor,'vox');
                stopTime = clock;
                
                clf
                %cellOri = plotSkeletons(cellArbors.arbors(i),seedSub);
                arbor = cellArbors.arbors(i);
                subplot(2,2,1)
                skelSholl = plotSholl(arbor,seedSub,[1 2]);
                skelSholls.YX = skelSholl;
                
                title('coronal')
                DOi(1) = skelSholl.DOi;
                subplot(2,2,2)
                skelSholl = plotSholl(arbor,seedSub,[1 3]);
                skelSholls.YZ = skelSholl;
                title('sagital')
                DOi(2) = skelSholl.DOi;
                subplot(2,2,3)
                skelSholl = plotSholl(arbor,seedSub,[2 3]);
                skelSholls.XZ = skelSholl;
                title('transverse')
                DOi(3) = skelSholl.DOi;
                
                
                
                
                
                subplot(4,2,6);
                bar(DOi)
                ylim([0 1])
                
                subplot(4,2,8);
                indPref = cellCon(:,i);
                prefColor = [indPref/max(indPref)];
                try
                    bar(indPref,'FaceColor', prefColor)
                end
                ylim([0 30])
                
                
                pause(.01)
                imageName = sprintf('%smorphSholl_%03.0f.png',shollConDirImage,cellTarg)
                dataName = sprintf('%smorphSholl_%03.0f.mat',shollConDirData,cellTarg)
                %saveas(gcf,fileName,'png')
                save(dataName,'skelSholls')
                
                
                
                
                foundSkel(i) = 1;
                
            else
                'cell Failed'
                
                foundSkel(i) = 0;
                clf
                pause(.01)
                imageName = sprintf('%smorphSholl_%03.0f.png',shollConDirImage,cellTarg)
                dataName = sprintf('%smorphSholl_%03.0f.mat',shollConDirData,cellTarg)
                saveas(gcf,fileName,'png')
                save(dataName,'skelSholls')
                
                
                
                
            end
            
        end
        
        
    end
    
    
    allDOi(i,:) = DOi;
end



%%

d = 1

subplot(1,1,1)

nanCell = sum(isnan(allDOi),2)>0;

useCell = (foundSkel>0)  & (isPref>2) & ~nanCell';
scatter(cellPref(useCell), allDOi(useCell,d))

xlim([-.2 1.2])
ylim([0 1])

clf
subplot(3,1,1)
hist(allDOi((useCell & (cellPref==0)),d))

subplot(3,1,2)
hist(allDOi((useCell & (cellPref>0) & (cellPref<1)),d))

subplot(3,1,3)
hist(allDOi((useCell & (cellPref==1)),d))










