
clear all

%MPN = GetMyDir

load('MPN.mat');
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])
useImage  = 1;

[allTraced tracedAxons tracedTCR] = getList_tracedCells;
cellList = tracedTCR;
[colIDs colProp] = getList_inSpine;
colProp = colProp/max(colProp)*4;

%% Skeletonize with shortest paths

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
seedList = [108 201 903 109 907];

useList = obI2cellList_seedInput(obI,seedList);
conTo = makeConTo(obI,seedList);
seedPref = seedPreferences(seedList,useList);

cellCon =  seedPref.cellPrefNoExclusion;
cellPref = seedPref.cellPrefNoExclusion(1,:)./sum(seedPref.cellPrefNoExclusion(1:2,:));
isPref = sum(seedPref.cellPrefNoExclusion(1:2,:));

%%Pick seed group
sharedSyn = seedPref.sharedSyn;
seedGroupAll = zeros(1,size(sharedSyn,2));
for i = 1:size(sharedSyn,2);
   seedCon = sharedSyn(1:2,i)>0; %% only group to first two seed cells
   if sum(seedCon) == 1
       seedGroupAll(i) = find(seedCon);
   end    
end

[A idx] = intersect(useList.postList,cellList);
seedGroup = seedGroupAll(idx);

%% Get cell images
if useImage
    Dim = 1;
    dsampI = .1;
    fsize = max(cat(1,dsObj(:).subs),[],1);
    anchorScale = [.0184 0.016 0.030];
    voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
    clear allYXV
    for i = 1:length(cellList)
        
        col = [1 1 1];
        Iraw = showCellSum(obI,dsObj,cellList(i),col,Dim,fsize);
        I = Iraw(:,:,1);
        I = imresize(I,dsampI,'bicubic');
        [y x] = find(I >0);
        v = I(I>0);
        peaks = find(v>median(v));
        my = median(y(peaks));
        mx = median(x(peaks));
        
        allYXV{i} = [y-my x - mx v];
    end
end



%%
c = 0;
clear recRes latRat cellArbors
dsamp = 1;
for i = 1:length(cellList)
    clf
    subplot(2,1,1)
    
    cellTarg = cellList(i);
    disp(sprintf('skeletonizing cell %d (%d of %d)',cellTarg,i,length(cellList)))
    
    skelFile = sprintf('%s%d.mat',TPNmat,cellTarg);
    viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    
    
    rawObjectSubs = getCellSubs(obI,dsObj,cellTarg);
    rawSeed   = ceil(getSeed(obI,cellTarg));
    
    if (size(rawObjectSubs,1)>100) & (~isempty(rawSeed))
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
            subplot(1,1,1)
            'cell passed'
            %             cellView = cellStruct.sideViews{1};
            %             skel = cellStruct.skel;
            lookRange = [5 25];
            [cellOri pc] = princompSkeletons(cellStruct,seedSub,lookRange);
            pause(.1)
            
            cellArbors(i).pc = pc;
            cellArbors(i).cellOri = cellOri;
            
            Score = pc.Score;
            yxView = round(Score(:,[1 2])/dsamp);
            uyxView = unique(yxView,'rows');
            %             scatter(uyxView(:,1),uyxView(:,2),'.')
            %             xlim([-50 50])
            %             ylim([-50 50])
            
            cellArbors(i).yxView = yxView;
            cellArbors(i).uyxView = uyxView;
            
            
            %             subplot(2,1,2)
            %             indPref = cellCon(:,i);
            %             prefColor = [indPref/max(indPref)];
            %             bar(indPref,'FaceColor', prefColor)
            %             ylim([0 30])
            %
            pause(.1)
            
            fileName = sprintf('%smorphOri_%03.0f.png',morphConDir,cellTarg)
            % saveas(gcf,fileName,'png')
            
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
            
            foundSkel(i) = 0;
            
        end
        
        
    end
    
end

if sum(foundSkel==0)
    disp('missing skeleton')
end

%%


cellPrefRaw = seedPref.sharedAx;
cellPref(1:2,:) = cellPrefRaw(1:2,:);
cellPref(3,:) = max(cellPrefRaw(3:end,:),[],1);
prefList = seedPref.cellList;
prefScale = 10;

%%
%%Pick property
clear cellProps otherProp
for i = 1:length(cellArbors)
    
    lats = cellArbors(i).cellOri.latent;
    %lats = cellArbors(i).cellOri.oriCount;
    prop = [lats(2)/lats(1) lats(3)/lats(2)];
    %prop = [(lats(1)-lats(2))/lats(1) (lats(1)-lats(3))/lats(1)];
    %prop = [lats(2) lats(3)];
    otherProp(i,1) = [lats(3)/lats(1)];
    %prop = cellArbors(i).cellOri.latRat;
    cellProps(i,:)  =prop;
    keepLats(i,:) = lats;
end

%cellProps = cellProps/max(keepLats(:,2));

%%Make histograms
histBin = [0:.05:1];
subplot(3,1,1)
histProp1 = histc(cellProps(:,1),histBin)
bar(histBin,histProp1,'barwidth',1)
xlim([0 1])
subplot(3,1,2)
histProp2 = histc(cellProps(:,2),histBin)
bar(histBin,histProp2,'barwidth',1)
xlim([0 1])

subplot(3,1,3)
histProp3 = histc(otherProp,histBin)
bar(histBin,histProp3,'barwidth',1)
xlim([0 1])

%% Analyze group 1 and 2
cellProps1 = cellProps(seedGroup==1,:);
mean1 = mean(cellProps1,1);
cellProps2 = cellProps(seedGroup==2,:);
mean2 = mean(cellProps2,1);
ranksum(cellProps1(:,1),cellProps2(:,1))



%%Bootstrap 1 and 2
realD2mean = sqrt((mean1(1)-mean2(1)).^2 + (mean2(2)-mean1(2)).^2);
allProps = [cellProps1; cellProps2];
N = size(allProps,1);
N1 = size(cellProps1,1); N2 = size(cellProps2,1);
randD2mean = zeros(10000,1);
for r = 1:10000
   newProp = allProps(randperm(N),:);
   newMean1 = mean(newProp(1:N1,:),1);
   newMean2 = mean(newProp(N1+1:end,:),1);
   randD2mean(r) = sqrt((newMean1(1)-newMean2(1)).^2 + (newMean1(2)-newMean2(2)).^2);
end

histBinD2mean = [0:.01:1];
histMCdist2mean = hist(randD2mean,histBinD2mean);
bar(histBinD2mean,histMCdist2mean,'barwidth',1)
hold on
scatter(realD2mean,.1,40,'r','filled')
hold off

mean(randD2mean)
realD2mean
P = sum(randD2mean>=realD2mean)/10000
sortRand = sort(randD2mean);
range95 = [sortRand(round(.025 * 10000)) sortRand(round((1-.025) * 10000))] 

%% Analyze group 1 and 2 - OtherProp
cellProps1 = otherProp(seedGroup==1);
mean1other = mean(cellProps1);
cellProps2 = otherProp(seedGroup==2);
mean2other = mean(cellProps2);
ranksum(cellProps1,cellProps2)



%%Bootstrap 1 and 2- Other Prop
realD2mean = abs(mean1other - mean2other)
allProps = [cellProps1; cellProps2];
N = size(allProps,1);
N1 = size(cellProps1,1); N2 = size(cellProps2,1);
randD2mean = zeros(10000,1);
for r = 1:10000
   newProp = allProps(randperm(N));
   newMean1 = mean(newProp(1:N1));
   newMean2 = mean(newProp(N1+1:end));
   randD2mean(r) = abs(newMean1 - newMean2);
end

histBinD2mean = [0:.01:1];
histMCdist2mean = hist(randD2mean,histBinD2mean);
bar(histBinD2mean,histMCdist2mean,'barwidth',1)
hold on
scatter(realD2mean,.1,40,'r','filled')
hold off

P = sum(randD2mean>=realD2mean)/10000
sortRand = sort(randD2mean);
mean(randD2mean)
mean(sortRand)
range95 = [sortRand(round(.025 * 10000)) sortRand(round((1-.025) * 10000))] 

%% Bar used for model
histBin = [0:.05:1];
subplot(3,1,1)
histProp1 = histc(cellProps(seedGroup>0,1),histBin)
bar(histBin,histProp1,'barwidth',1)
xlim([0 1])
subplot(3,1,2)
histProp2 = histc(cellProps(seedGroup>0,2),histBin)
bar(histBin,histProp2,'barwidth',1)
xlim([0 1])

subplot(3,1,3)
histProp3 = histc(otherProp(seedGroup>0),histBin)
bar(histBin,histProp3,'barwidth',1)
xlim([0 1])

%% Kmeans grouping
% for r = 1:10
%     k = 4;
%     [IDX,C,sumd,D] = kmeans(cellProps,k);
%     groupCol = {'r','g','b','c','m','y'}
%     
%     clf,
%     for g = 1:k
%         scatter(cellProps(IDX==g,1),cellProps(IDX == g,2),groupCol{g},'filled');
%         hold on
%         
%     end
%     hold off
%     pause(.05)
% end

spread = 2000;
buf = 100;
if useImage
    maxVals = max(abs(cat(1,allYXV{:})),[],1);
    buf = max(ceil(maxVals(1:2)));
end
field = zeros(spread+ buf*2,spread+ buf*2,3);
clf
hold on
for i = 1:length(cellArbors)
    cellProp = cellProps(i,:);%cellArbors(i).cellOri.latRat;
    
    cellShift = round(cellProp * spread);
    cellShift(cellShift<1) = 1;
    cellShift(cellShift>spread) = spread;
    
    
    if useImage
        clear yxPos
        yxv = allYXV{i};
        yxPos(:,1) = yxv(:,1) + cellShift(:,1)+buf;
        yxPos(:,2) = cellShift(:,2) + yxv(:,2)+buf;
        yxVal = yxv(:,3)/max(yxv(:,3));
        
    else
        %cellProp = abs(cellArbors(i).cellOri.oriRat);
        %
        %
        %
        %     yxPos = yxView*0;
        %     yxPos(:,1) = cellShift(:,1) + yxView(:,1);
        %     yxPos(:,2) = cellShift(:,2) + yxView(:,2);
        %     scatter(yxPos(:,1),yxPos(:,2),16,'.')
        
        
        yxView = cellArbors(i).uyxView;
        
        yxView = cellArbors(i).yxView;
        yxPos = yxView*0;
        yxPos(:,1) = cellShift(:,1) + yxView(:,1) + buf;
        yxPos(:,2) = cellShift(:,2) + yxView(:,2) + buf;
        yxVal = ones(size(yxPos,1));
    end
    
    yxPos = ceil(yxPos);
    
    if 0 %plot by seed connectivity
    prefTarg = find(prefList == cellList(i));
    if ~isempty(prefTarg)
        addVals = cellPref(:,prefTarg)*prefScale;
    else
        addVals = [1; 1; 1];
    end
    else
       prefTarg = find(colIDs == cellList(i));
       
       if isempty(prefTarg)
           addVals = [.3; .3; .3]
       else
           addVals = [colProp(prefTarg) 0 1-colProp(prefTarg)];
       end
        
        
    end
    addVals(addVals<0) = 0;
    
    addVals = addVals/max(addVals);
    for p = 1:size(yxPos,1)
        for c = 1:3
            addVal = addVals(c);
            field(yxPos(p,2),yxPos(p,1),c) =  field(yxPos(p,2),yxPos(p,1),c) + addVal * yxVal(p);
        end
    end
end
%field = field(buf+1:buf+spread, buf+1:buf+spread,:);
ylim([0 spread])
xlim([0 spread])
hold off

image(uint8(field*600))

% imwrite(uint8(field*10),'D:\LGNs1\Analysis\dendriticMorphology\out1\sharedAxVsPrinComp2.png')
%%%
%%
imshow(uint8(field*1200))
hold on

for i = 1:length(cellArbors)
    
    yxView = cellArbors(i).uyxView;
    cellProp = cellProps(i,:); %cellArbors(i).cellOri.latRat;
    %cellProp = abs(cellArbors(i).cellOri.oriRat);
    
    cellShift = round(cellProp * spread)+buf+spread/100;
    cellShift(cellShift<1) = 1;
    cellShift(cellShift>spread) = spread;
    
    cellText = sprintf('%.2f',otherProp(i));
    %text(cellShift(1),cellShift(2),num2str(cellList(i)),'color','w','fontsize',8)
    %text(cellShift(1),cellShift(2),num2str(IDX(i)),'color','w','fontsize',8)
   % text(cellShift(1),cellShift(2),cellText,'color','w','fontsize',8)

end


posMean1 = round(mean1 * spread)+buf+spread/100;
posMean2 = round(mean2 * spread)+buf+spread/100;
text(posMean1(1),posMean1(2),'*','color','r','fontsize',50)
text(posMean2(1),posMean2(2),'*','color','g','fontsize',50)


plot([buf spread+buf],[buf buf],'w')
plot([buf buf],[buf spread+buf],'w')
text(buf/2,buf/2,num2str(0),'color','w')
text(buf/2,buf+spread,'third','color','w')
text(buf+spread,buf/2,'second','color','w')

%{

epsName = 'D:\LGNs1\Analysis\dendriticMorphology\out5\labeledLatRat.eps'
print(gcf, epsName, '-depsc2','-painters')
pngName = 'D:\LGNs1\Analysis\dendriticMorphology\out4\labeledLatRat.png'
print(gcf,pngName,'-dpng','-r1024','-opengl','-noui')
imwrite(field,'D:\LGNs1\Analysis\dendriticMorphology\out5\field.png')
%}
hold off
%%










