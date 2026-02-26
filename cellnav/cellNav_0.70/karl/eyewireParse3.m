%% Setup
loadBool=1;
%set the source directories
jsonDirBPC='G:\Data\eyewire\bpc\';
jsonDirRGC='G:\Data\eyewire\rgc\';
%get the file lists
jsonDirStructBPC=dir(jsonDirBPC);
jsonDirStructRGC=dir(jsonDirRGC);
%trim to just json files
jsonFileListBPC=string([]);
jsonFileListRGC=string([]);
for k=1:length(jsonDirStructBPC)
    curFileName=string(jsonDirStructBPC(k).name);
    if endsWith(curFileName,'json')
        jsonFileListBPC=[jsonFileListBPC; curFileName];
    end
end

for l=1:length(jsonDirStructRGC)
    curFileName=string(jsonDirStructRGC(l).name);
    if endsWith(curFileName,'json')
        jsonFileListRGC=[jsonFileListRGC; curFileName];
    end
end

%% Put everything into a structure
%load in bpc by type

if loadBool==1
    %make a struct to hold all the bpc data
    global cellDat;
    cellDat=struct;
    global ewIDList;
    ewIDList=[];
    cellIt=1;
    %load bpc id, type, and stratDat
    for k=1:length(jsonFileListBPC)
        curFileName=[jsonDirBPC+jsonFileListBPC(k)];
        rawTxt=fileread(curFileName);
        curJ=jsondecode(rawTxt);
        for curCell=1:length(curJ)
            curCellDat=curJ(curCell);
            cellDat(cellIt).id=curCellDat.id;
            ewIDList=[ewIDList curCellDat.id];
            cellDat(cellIt).type=curCellDat.type;
            cellDat(cellIt).strat=curCellDat.stratification;
            cellIt=cellIt+1;
        end
    end
    
    global rgcDat;
    rgcDat = struct;
    global ewRGCList;
    ewRGCList=[];
    cellIt=1;
    for k=1:length(jsonFileListRGC)
        curFileName=[jsonDirRGC+jsonFileListRGC(k)];
        rawTxt=fileread(curFileName);
        curJ=jsondecode(rawTxt);
        for curCell=1:length(curJ)
            curCellDat=curJ(curCell);
            rgcDat(cellIt).id=curCellDat.id;
            ewRGCList=[ewRGCList curCellDat.id];
            rgcDat(cellIt).type=curCellDat.type;
            rgcDat(cellIt).strat=curCellDat.stratification;
            cellIt=cellIt+1;
        end
    end
    
    %% Josh Setup
    global fvDir;
    fvDir='G:\Data\eyewire\Analysis\cellNavDat\';
    fvDir='Y:\karlsRetina\Analysis\cellNav\Karl_210121\DF2020\Volumes\DF2020_2\Analysis\fvLibrary\';
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\fvLibrary\';
    plotDest=[fvDir '\dat\'];
    if ~exist(plotDest,'dir')
        mkdir(plotDest)
    end
    %tis = load([fvDir 'tis.mat']);
    
    %These are necessary for the IPLdepth function to run
    if exist([fvDir 'ref_gcl nucEdge.mat'],'file')
        ipl_bord_GCL = load([fvDir 'ref_gcl nucEdge.mat']);
        ipl_bord_INL = load([fvDir 'ref_inl nucEdge.mat']);
        %fit some planes to the borders of the IPL and INL
        GCLbord=ipl_bord_GCL.fv.vertices(:,:);
        INLbord=ipl_bord_INL.fv.vertices(:,:);
        
    else
        GCLbord = [0 0 0; 100 0 0; 0 100 0; 100 100 0];
        INLbord = [0 0 100; 100 0 100; 0 100 100; 100 100 100];
    end
    
    newMethod=1;
    if newMethod==0
        %convert border vertices to pt cloud
        GCLptCloud=pointCloud(GCLbord);
        INLptCloud=pointCloud(INLbord);%fit plane to the pt clouds (result params = ax+by+cz+d=0)
        GCLplane=pcfitplane(GCLptCloud,1,'MaxNumTrials',50000);
        INLplane=pcfitplane(INLptCloud,1,'MaxNumTrials',50000);
    elseif newMethod==1
        Locs={GCLbord;INLbord}; %These are in z,x,y I think.
        for i=1:2
            P=Locs{i};
            B(:,i) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
        end
        GCLplane=struct();
        INLplane=struct();
        GCLplane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
        INLplane.Parameters=[-1 B(2,2) B(1,2) B(3,2)];
    end
    
    %get the cids of all bpcs
    load([fvDir 'tis.mat']);
    global bpcCidList
    bpcCidList=tis.cids;
    
    %% Loop through all of the bipolar cells and get their depths within the IPL
    
    %save individual histograms of all bpc cells
    saveJPGBool=0;
    
    %create empty struct for results
    global bpcStats;
    bpcStats=struct;
    
    for i=1:length(bpcCidList)
        %get the cid
        bpcStats(i).cid=bpcCidList(i);
        
        matString=sprintf('%i.mat',bpcStats(i).cid);
        if ~exist([fvDir matString],'file')
            bpcStats(i).isDat = 0;
        else
            matDat=load([fvDir matString]);
            if isempty(matDat.fv.vertices)
                bpcStats(i).isDat = 0;
            else
                bpcStats(i).isDat = 1;
                %get the depth in the IPL
                bpcStats(i).zDepth=getCellZdepth(bpcCidList(i),GCLplane,INLplane);
                %go ahead and get the old percentage for comparison
                %bpcStats(i).zDepthOld=getCellZdepthOld(bpcCidList(i),ipl_bord_GCL,ipl_bord_INL);
                %get area modeled with ellipse
                bpcStats(i).area=getCellArborArea(bpcCidList(i));
                %get average span in x / y
                bpcStats(i).spanMicron=getCellArborSpan(bpcCidList(i));
                %get the average vertex distance from centroid
                [bpcStats(i).densityAvg,bpcStats(i).densityStd]=getCellArborDensity(bpcCidList(i));
                %JOSH; HERE IS THE HISTOGRAM FUNCTION
                %get the data for plotting the histograms
                bpcStats(i).histDat=getDepthDist(bpcCidList(i),GCLplane,INLplane,saveJPGBool,plotDest);
            end
        end
    end
    
    %% get histogram data for doing the comparison things.
    fetchAll=1;
    if fetchAll==1
        superCidList=[bpcCidList cellDat(:).id];
        for cellIt=1:length(superCidList)
            curCid=superCidList(cellIt);
            if sum(ismember(bpcCidList,curCid))
                %get the histogram for our cell
                curHistDatRaw=bpcStats(bpcCidList==curCid).histDat;
                if ~isempty(curHistDatRaw)
                    curHistDatRaw(:,2)=curHistDatRaw(:,2);
                end
            elseif sum(ismember(ewIDList,curCid))
                %get the histogram for eyewire cell
                curHistDatRaw=cellDat(ewIDList==curCid).strat;
                curHistDatRaw(:,2)=curHistDatRaw(:,2);
            end
            if ~isempty(curHistDatRaw) && sum(curHistDatRaw(:,2))>2000
                %normalize the histogram
                normHistDat=curHistDatRaw(:,2)./(sum(curHistDatRaw(:,2)*1));
                %plot the histogram
                %plot(curHistDatRaw(:,1),normHistDat);
            end
        end
    end
    
    %% Just the data; no figs BROKEN
    
    runDat=0;
    if runDat==1
        allTheCids=[bpcStats(:).cid cellDat(:).id];
        %h=jumboHisto(allTheCids,0,string("asdf"));
        
    end
    
    
    %% GraphDat creation
    graphDat=struct;
    graphDat(1).bcSubtype=3;
    graphDat(1).bcEW='bc3a';
    graphDat(2).bcSubtype=4;
    graphDat(2).bcEW='bc3b';
    graphDat(3).bcSubtype=5;
    graphDat(3).bcEW='bc4';
    graphDat(4).bcSubtype=6;
    graphDat(4).bcEW='bc5i';
    graphDat(5).bcSubtype=7;
    graphDat(5).bcEW='bc5o';
    graphDat(6).bcSubtype=8;
    graphDat(6).bcEW='bc5t';
    
end
%% Get the cells of a certain type
bipolDatStruct=struct;
runNorm=1;
if runNorm==1
    for i=1:length(graphDat)
        curCidList=[];
        for curCellIt=1:length(tis.cells.cids)
            if tis.cells.type.typeID(curCellIt)==7 && tis.cells.type.subTypeID(curCellIt)==graphDat(i).bcSubtype
                curCidList=[curCidList tis.cells.cids(curCellIt)];
            end
        end
        for curCellIt=1:length(cellDat)
            if string(getfield(cellDat(curCellIt),'type'))==graphDat(i).bcEW
                curCidList=[curCidList double(getfield(cellDat(curCellIt),'id'))];
            end
        end
        
        curTitle=graphDat(i).bcEW;
        
        bipolDatStruct(i).dat=jumboHisto('bpc',curCidList,0,string(curTitle));
        bipolDatStruct(i).name=graphDat(i).bcEW;
        %h=bipolDatStruct.dat;
    end
end

%% test out the getStats function
getTheStats=1;
if getTheStats==1
    allStatDat=struct;
    for i=1:length(graphDat)
        curDat=bipolDatStruct(i).dat;
        allStatDat(i).dat=getStats(curDat);
    end
end


%% revised Fitting
figure;
hold on
testBool=1;
if testBool==1
    for bpcType=1:length(bipolDatStruct)
        curTypeDat=bipolDatStruct(bpcType).dat;
        curTypeMeanML=[];
        curTypeMeanEW=[];
        for i=1:length(curTypeDat)
            if curTypeDat(i).source==1
                %ew bpc
                curCellHistDat=curTypeDat(i).histDat;
                curCellArray=[];
                for rowIt=1:length(curCellHistDat)
                    if curCellHistDat(rowIt,2)>0
                        curCellArray=[curCellArray curCellHistDat(rowIt,1)*curCellHistDat(rowIt,2)];
                    end
                end
                %curCellMean=mean(curCellArray>0);
                hMax=max(curCellHistDat(:,2));
                CHDC=curCellHistDat;
                CHDC(CHDC(:,2)<(hMax*.3),2)=0;
                curAvg=sum(CHDC(:,1).*CHDC(:,2))/sum(CHDC(:,2));
                curTypeMeanEW=[curTypeMeanEW curAvg];
            elseif curTypeDat(i).source==2
                %morg bpc
                curCellHistDat=curTypeDat(i).histDat;
                curCellArray=[];
                for rowIt=1:length(curCellHistDat)
                    if curCellHistDat(rowIt,2)>0
                        curCellArray=[curCellArray curCellHistDat(rowIt,1)*curCellHistDat(rowIt,2)];
                    end
                end
                %curCellMean=mean(curCellArray>0);
                hMax=max(curCellHistDat(:,2));
                CHDC=curCellHistDat;
                CHDC(CHDC(:,2)<(hMax*.3),2)=0;
                curAvg=sum(CHDC(:,1).*CHDC(:,2))/sum(CHDC(:,2));
                curTypeMeanML=[curTypeMeanML curAvg];
            end
            bipolDatStruct(bpcType).mlMeans=curTypeMeanML;
            bipolDatStruct(bpcType).ewMeans=curTypeMeanEW;
            
            %total arrays for each type
            %get means and SDs for each
            %put results into results mat
            %graph results
        end
    end
end

%% scatter cell means
figure;
hold on;
for bpcType=1:length(bipolDatStruct)
    scatter(rand(length(bipolDatStruct(bpcType).ewMeans),1)+bpcType*2, bipolDatStruct(bpcType).ewMeans, 'r');
    scatter(rand(length(bipolDatStruct(bpcType).mlMeans),1)+bpcType*2, bipolDatStruct(bpcType).mlMeans, 'g');
    
end
names = {'bpc3a'; 'bpc3b'; 'bpc4'; 'bpc5i'; 'bpc5o'; 'bpc5t'};
set(gca,'xtick',[1:6],'xticklabel',names)

%% scatter type means
figure;
hold on;
cols=['r','g','b','m','c','k'];
for bpcType=1:length(bipolDatStruct)
    curTypeEwMean=mean(bipolDatStruct(bpcType).ewMeans);
    curTypeCol=cols(bpcType);
    for mlBpcIt=1:length(bipolDatStruct(bpcType).mlMeans)
        curMlBpcMean=bipolDatStruct(bpcType).mlMeans(mlBpcIt);
        curOffset=curMlBpcMean-curTypeEwMean;
        bipolDatStruct(bpcType).offSet(mlBpcIt)=curOffset;
        scatter(curMlBpcMean,curOffset,curTypeCol);
    end
end



%% scatter type means
scaleSearch=zeros(20,28);
targVal = []; realVal = [];
for bpcType=1:length(bipolDatStruct)
    curTypeEwMean=mean(bipolDatStruct(bpcType).ewMeans);
    curTypeCol=cols(bpcType);
    for mlBpcIt=1:length(bipolDatStruct(bpcType).mlMeans)
        curMlBpcMean=bipolDatStruct(bpcType).mlMeans(mlBpcIt);
        curOffset=curMlBpcMean-curTypeEwMean;
        bipolDatStruct(bpcType).offSet(mlBpcIt)=curOffset;
        targVal = [targVal curTypeEwMean];
        realVal = [realVal curMlBpcMean];
    end
end

realErr = targVal - realVal;
totRealErr = sum(abs(realErr));
%%
%hold off
testScale = [.96:.0001:.98];
testAdd = [-.03:.0001:-.02];
totTestErr = zeros(length(testScale),length(testAdd));
for s = 1:length(testScale)
    for a = 1:length(testAdd)
        
        testVal = realVal * testScale(s) + testAdd(a);
        testErr = targVal - testVal;
        totTestErr(s,a) = sum(abs(testErr));
    end
end
colormap jet(256)
imErr = totTestErr;
imErr = imErr - min(imErr(:));
imErr = imErr/max(imErr(:));
imErr = 1-imErr;
%image(imErr * 256),pause(.01)
[bestS bestA] = find(totTestErr == min(totTestErr(:)));
bestScale = testScale(bestS)
bestAdd = testAdd(bestA)

%% real math
n=28;
x=realVal;
y=targVal;
m=(n*sum(x.*y)-(sum(x)*sum(y)))/(n*sum(x.^2)-(sum(x)^2))
b=(sum(y)-m*sum(x))/28

%% find the fit of eyewire to us
testBool=0;
if testBool==1
    outputStruct=struct;
    comparisonDat=zeros(length(graphDat),4);
    outputMatDatStruct=struct;
    figure();
    hold on
    for i=1:length(graphDat)
        %outputMatDatStruct=struct;
        curCidList=[];
        ewMorg=[];
        for curCellIt=1:length(tis.cells.cids)
            if tis.cells.type.typeID(curCellIt)==7 && tis.cells.type.subTypeID(curCellIt)==graphDat(i).bcSubtype
                curCidList=[curCidList tis.cells.cids(curCellIt)];
                ewMorg=[ewMorg 2];
            end
        end
        for curCellIt=1:length(cellDat)
            if string(getfield(cellDat(curCellIt),'type'))==graphDat(i).bcEW
                curCidList=[curCidList double(getfield(cellDat(curCellIt),'id'))];
                ewMorg=[ewMorg 1];
            end
        end
        outputMatDatStruct(i).cid=curCidList;
        %[outputStruct.cid];
        outputMatDat=zeros(length(curCidList),3);
        for cellIt=1:length(curCidList)
            curCid=curCidList(cellIt);
            dataSet=ewMorg(cellIt);
            targField=find(cidListHistDat==curCid);
            outputMatDat(cellIt,1)=curCid;
            outputMatDat(cellIt,2)=outputStruct(targField).histMean;
            outputMatDat(cellIt,3)=outputStruct(targField).histSD;
        end
        comparisonDat(i,1)=mean(outputMatDat(ewMorg==1,2));
        comparisonDat(i,2)=mean(outputMatDat(ewMorg==2,2));
        comparisonDat(i,3)=mean(outputMatDat(ewMorg==1,3));
        comparisonDat(i,4)=mean(outputMatDat(ewMorg==2,3));
        scatter(i,comparisonDat(i,1),'r');
        errorbar(i,comparisonDat(i,1),comparisonDat(i,3),'r');
        scatter(i,comparisonDat(i,2),'b');
        errorbar(i,comparisonDat(i,2),comparisonDat(i,4),'b');
        outputMatDatStruct(i).data=outputMatDat;
    end
end




%% Fitting room

%Getting the RGC histo overlaps

%Everything is already loaded into the bpcStats, so we just need to get our
%RGC list to iterate through.
cidListRGC=tis.cells.cids(tis.cells.type.typeID==1);

RGCgraphDat=struct;
rgcTypeList=string('1ni');
for i=1:length(rgcDat)
    curType=string(rgcDat(i).type);
    if ~ismember(curType,rgcTypeList)
        rgcTypeList=[rgcTypeList curType];
    end
end
for j=1:length(rgcTypeList)
    curType=rgcTypeList(j);
    RGCgraphDat(j).bcEW=curType;
end

RGCDatStruct=struct;
runNorm=1;
if runNorm==1
    for i=1:length(RGCgraphDat)
        curCidList=[];
        for curCellIt=1:length(tis.cells.cids)
            if tis.cells.type.typeID(curCellIt)==999 %&& tis.cells.type.subTypeID(curCellIt)==graphDat(i).bcSubtype
                curCidList=[curCidList tis.cells.cids(curCellIt)];
            end
        end
        for curCellIt=1:length(rgcDat)
            if string(getfield(rgcDat(curCellIt),'type'))==RGCgraphDat(i).bcEW
                curCidList=[curCidList double(getfield(rgcDat(curCellIt),'id'))];
            end
        end
        
        curTitle=RGCgraphDat(i).bcEW;
        
        RGCDatStruct(i).dat=jumboHisto('rgc',curCidList,0,string(curTitle));
        RGCDatStruct(i).name=RGCgraphDat(i).bcEW;
        %h=bipolDatStruct.dat;
    end
end

RGCDatStruct=getMeanHistos(RGCDatStruct,0);

%Go through our dataset and get all of the rgcs and make ADJUSTED
%histograms for them
m=0.97;
b=-0.02;
sizeCut=20000;
crossCorrResults=zeros(length(RGCDatStruct),length(cidListRGC));

for i=1:length(cidListRGC)
    curCid=cidListRGC(i);
    curHistDat=bpcStats(bpcCidList==curCid).histDat;
    %make sure that the cell isn't tiny
    if length(curHistDat)>0
        if sum(curHistDat(:,2))>sizeCut
            newZ=(curHistDat(:,1).*m)+b;
            percHist=curHistDat(:,2)./sum(curHistDat(:,2))*100;
            h=jumboHisto('rgc',curCid,0,string(curCid));
            h.histDat(:,2)=h.histDat(:,2).*100;
            for k=1:length(RGCDatStruct)
                curRGCTypeDat=RGCDatStruct(k).meanHist;
                curRGCTypeName=RGCDatStruct(k).name;
                offSet=abs(curRGCTypeDat(:,1)-h.histDat(1,1));
                lag=find(offSet==min(offSet));
                emptySet=curRGCTypeDat.*[1 0];
                emptySet(lag:lag+length(h.histDat(:,1))-1,2)=h.histDat(:,2);
                debugYes=0;
                if debugYes==1
                    figure(); hold on
                    scatter(h.histDat(:,1),h.histDat(:,2))
                    scatter(emptySet(:,1),emptySet(:,2))
                    scatter(curRGCTypeDat(:,1),curRGCTypeDat(:,2))
                end
                %remove everything above 0.95 for getting rid of soma
                somaCut=0.95;
                emptySet(emptySet(:,1)>somaCut,2)=0;
                curCorr=corr(emptySet(:,2),curRGCTypeDat(:,2));
                crossCorrResults(k,i)=curCorr;
            end
        end
    end
end

%%

figure(); hold on
cmap=jet(256);
image(crossCorrResults.*256);
colormap(cmap);

names = string(cidListRGC);
set(gca,'xtick',[1:length(names)],'xticklabel',names)
xtickangle(90)

set(gca,'ytick',[1:length(rgcTypeList)],'yticklabel',rgcTypeList)
%%
compPix=[40 41 50 65];
compPix=[9 24 16];
refPix=[19];
figure();
hold on
for i=1:length(compPix)
    curCid=cidListRGC(compPix(i));
    %title(string(curCid));
    curHistDat=bpcStats(bpcCidList==curCid).histDat;
    newZ=(curHistDat(:,1).*m)+b;
    somaCut=0.9;
    curHistDat(newZ>somaCut,2)=0;
    percHist=curHistDat(:,2)./sum(curHistDat(:,2))*100;
    percHistSmooth=smoothdata(percHist,1,'movmean',3);
    %plot(curHistDat(:,1),percHist,'r');
    plot(curHistDat(:,1),percHistSmooth);
    if i>1
        names=[names string(curCid)];
    else
        names=[string(curCid)];
    end
end
for i=1:length(refPix)
    curRGCTypeDat=RGCDatStruct(refPix(i)).meanHist;
    names=[names RGCDatStruct(refPix(i)).name];
    plot(curRGCTypeDat(:,1),curRGCTypeDat(:,2));
end

legend(names);


%% Function block

%get the average histos for each type in a jumbo output struct.
function outputStruct=getMeanHistos(inputStruct,plotBool)
outputStruct=inputStruct;
for i=1:length(inputStruct)
    curTypeHistDat=inputStruct(i).dat;
    totalHist=curTypeHistDat(1).histDat;
    for j=1:length(curTypeHistDat)
        curCellHist=curTypeHistDat(j).histDat;
        curCellHistDatPerc=curCellHist(:,2)./sum(curCellHist(:,2))*100;
        outputHist=curCellHist;
        outputHist(:,2)=curCellHistDatPerc;
        totalHist(:,2)=totalHist(:,2)+outputHist(:,2);
    end
    totalHist(:,2)=totalHist(:,2)./length(curTypeHistDat);
    outputStruct(i).meanHist=totalHist;
    if plotBool
        figure();
        hold on
        plot(totalHist(:,1),totalHist(:,2),'r');
        title(string(inputStruct(i).name));
    end
end
end

%get the stats for the histograms
function outputStruct=getStats(inputStruct)
outputStruct=inputStruct;
for curIt=1:length(inputStruct)
    curHistDat=inputStruct(curIt).histDat;
    curVertMat=[];
    for curBinIt=1:length(curHistDat)
        curBinZ=curHistDat(curBinIt,1);
        curBinVert=curHistDat(curBinIt,2)*10000;
        curVertMat=[curVertMat repmat(curBinZ,1,round(curBinVert))];
    end
    outputStruct(curIt).vertMat=curVertMat;
    outputStruct(curIt).histMean=mean(curVertMat);
    outputStruct(curIt).histSD=std(curVertMat);
end
%return outputStruct;
end

%make a histogram with a bunch of cells on it.
function histerical=jumboHisto(type,cidList,plotBool,titleString)
global bpcCidList;
global cellDat;
global bpcStats;
global ewIDList;
global ewRGCList;
global rgcDat;
if type=='bpc'
    useDat=cellDat;
    useList=ewIDList;
elseif type=='rgc'
    useDat=rgcDat;
    useList=ewRGCList;
end
histerical=struct;
if plotBool==1
    figure();
    hold on;
end
for cellIt=1:length(cidList)
    curCid=cidList(cellIt)
    histerical(cellIt).cid=curCid;
    if sum(ismember(bpcCidList,curCid))
        %get the histogram for our cell
        histerical(cellIt).source=2;
        curHistDatRaw=bpcStats(bpcCidList==curCid).histDat;
        if ~isempty(curHistDatRaw)
            curHistDatRaw(:,2)=curHistDatRaw(:,2);
            curHistDatRaw(:,2)=smoothdata(curHistDatRaw(:,2));
        end
        curHistSource=2;
    elseif sum(ismember(useList,curCid))
        histerical(cellIt).source=1;
        %get the histogram for eyewire cell
        curHistDatRaw=useDat(useList==curCid).strat;
        curHistDatRaw(:,2)=curHistDatRaw(:,2)*10;
        curHistSource=1;
    end
    if ~isempty(curHistDatRaw) && sum(curHistDatRaw(:,2))>2000
        %normalize the hiram
        normHistDat=curHistDatRaw(:,2)./(sum(curHistDatRaw(:,2)*1));
        %plot the histogram
        outputDat=[curHistDatRaw(:,1) normHistDat];
        
        histerical(cellIt).histDat=outputDat;
        if plotBool==1
            if curHistSource==1
                curCol='r';
            else
                curCol='g';
            end
            plot(curHistDatRaw(:,1),normHistDat,curCol);
        end
    end
end
if plotBool==1
    legend(string(cidList))
    ylim([0 .04]);
    xlim([-0.25 1.25]);
    title(titleString);
end
end

function plotDat=getDepthDist(cid,GCLplane,INLplane,plotBool,plotLoc)
global fvDir
%find the gcl and inl boundaries for percentage calcs
location=getMedianLoc(cid);
[gclz,inlz,zperc]=getIPLdepth(location(1),location(2),location(3),GCLplane,INLplane);
iplDepth = abs(abs(inlz)-abs(gclz));
%bin width for the output histogram atam
binWidth=iplDepth/515;
%load the mat file for the cell
matString=sprintf('%i.mat',cid);
matDat=load([fvDir matString]);
%get the vertices
coords=matDat.fv.vertices(:,:);
%make the histogram so we can get the data
figure();
hist=histogram(coords(:,1),'BinWidth',binWidth);

%create the histogram data for later
plotDatX=hist.BinEdges;
%bin edges are longer than bin counts so remove the last one
plotDatX=plotDatX(1:end-1);
%make the array the center of each bin
plotDatX=plotDatX+(binWidth/2);
%convert to depth percentages using the gcl and inl depths
%plotDatXpercs=abs(plotDatX-abs(inlz))/abs(abs(inlz)-abs(gclz));
plotDatXpercs=-(plotDatX-inlz)/abs(abs(inlz)-abs(gclz));

%vertex counts at each z location
plotDatY=hist.BinCounts;
%OPTIONAL: convert to a percentage of total arbor verts
plotDatYpercs=plotDatY/sum(plotDatY);
close;
%create the exported array (change to plotDatYpercs if you want).
plotDat=[plotDatXpercs',plotDatY'];
%this loop saves jpgs of each histogram
if plotBool==1
    plotFilename=sprintf('%i.jpg',cid);
    plotFilePath=[plotLoc plotFilename];
    figure();
    plot(plotDatXpercs,plotDatYpercs);
    xlim([0 1]);
    camroll(-90);
    title(plotFilename);
    saveas(gcf,plotFilePath);
    close;
end

end

%get the mean distance and standard deviation of vertex dist from centroid
function [cellArborDensityAvg,cellArborDensityStd]=getCellArborDensity(cid)
global fvDir
matString=sprintf('%i.mat',cid);
matDat=load([fvDir matString]);
coords=matDat.fv.vertices(:,:);
cellMedianLoc=median(coords);
coordsTrimmed=coords(abs(coords(:,1)-cellMedianLoc(1))<std(coords(:,1)),:);
distances=sqrt((coordsTrimmed(:,2)-cellMedianLoc(2)).^2+(coordsTrimmed(:,3)-cellMedianLoc(3)).^2);
cellArborDensityAvg=mean(distances);
cellArborDensityStd=std(distances);
end

%get the ranges of the arbor vertices, convert to microns
function cellArborSpanAvg=getCellArborSpan(cid)
global fvDir
matString=sprintf('%i.mat',cid);
matDat=load([fvDir matString]);
coords=matDat.fv.vertices(:,:);
ranges=range(coords);
cellArborSpanAvg=mean([ranges(2) ranges(3)])/10.2;
end

%get the area as modeled by ellipse
function cellArborArea=getCellArborArea(cid)
global fvDir
matString=sprintf('%i.mat',cid);
matDat=load([fvDir matString]);
coords=matDat.fv.vertices(:,:);
ranges=range(coords);
cellArborArea=pi*ranges(2)*ranges(3);
end

%get the depth from the plane fits
function zDepth=getCellZdepth(cid,GCLplane,INLplane)
location=getMedianLoc(cid);
[gclz,inlz,zDepth]=getIPLdepth(location(1),location(2),location(3),GCLplane,INLplane);
end

%deprecated, just here for sanity checking
function zDepth=getCellZdepthOld(cid,ipl_bord_GCL,ipl_bord_INL)
location=getMedianLoc(cid);
groupSize=20;
zDepth=getIPLdepthOld(location(1),location(2),location(3),groupSize,ipl_bord_GCL,ipl_bord_INL);
end

%get the centroid of a cid
function cellMedianLoc = getMedianLoc(cid)
global fvDir
matString=sprintf('%i.mat',cid);
matDat=load([fvDir matString]);
cellMedianLoc=median(matDat.fv.vertices(:,:));
end

%get the depth of a point in the IPL from the fitted planes
function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCLplane,INLplane)
zGCL=(-x*GCLplane.Parameters(2)-y*GCLplane.Parameters(3)-GCLplane.Parameters(4))/GCLplane.Parameters(1);
zINL=(-x*INLplane.Parameters(2)-y*INLplane.Parameters(3)-INLplane.Parameters(4))/INLplane.Parameters(1);

IPLdepth=abs(z-zINL)/abs(zINL-zGCL);

end

%old method using the points in the border region to calculate.
%deprecated. just here for sanity checks.
function IPLdepth = getIPLdepthOld(z,x,y,groupSize,ipl_bord_GCL,ipl_bord_INL)
GCLbord=ipl_bord_GCL.fv.vertices(:,:);
INLbord=ipl_bord_INL.fv.vertices(:,:);
%get the distances
distancesGCL=sum((GCLbord(:,2:3)-[x,y]).^2,2);
distancesINL=sum((INLbord(:,2:3)-[x,y]).^2,2);
%get the sorted distance idxs
[distGCL,GCLidx]=sort(distancesGCL);
[distINL,INLidx]=sort(distancesINL);
%get mean z locatin of points in border
zGCL=mean(GCLbord(GCLidx(1:groupSize),1));
zINL=mean(INLbord(INLidx(1:groupSize),1));
%calculate percentage
IPLdepth=(zINL-z)/(zINL-zGCL);
end

