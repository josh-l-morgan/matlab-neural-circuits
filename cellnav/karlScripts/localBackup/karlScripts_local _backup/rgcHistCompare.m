%rgcHistCompare

%use eyewireParse3 to get the rgcDat structure before starting
jsonDirBPC='Y:\karlsRetina\eyewire_data\bpc\';
jsonDirRGC='Y:\karlsRetina\eyewire_data\rgc\';
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

%%
w3List=[17012;17035;17095;17098;17138;20037;26025;26039;26054;26085;26113;26136;26154;26177;20120;20153;20182;20212;20258];
cellStr=struct2cell(rgcDat);
typeList=cellStr(2,1,:);
typeList=squeeze(typeList);
typeStr="4ow";
idList=find(contains(typeList,typeStr));


w3IDs=find(ismember([rgcDat.id],w3List));
typeStratRaw=[rgcDat(idList).strat];
typeStratRaw=flip(typeStratRaw,1);
stratRaw=[rgcDat(find(ismember([rgcDat.id],w3List))).strat];
stratRaw=flip(stratRaw,1);
typeStratDat=typeStratRaw(:,[2:2:end]);
stratDat=stratRaw(:,[2:2:end]);

%% out of curiosity, how does the vertex distribution compare to the volume?
fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';

%% test out some plots
w3sum=sum(stratDat,2);
typeSum=sum(typeStratDat,2);

w3Interp=interp1(stratRaw(:,1),w3sum,[-0.2:0.01:1.19]);
typeInterp=interp1(typeStratRaw(:,1),typeSum,[-0.2:0.01:1.19]);

tsTest=timeseries(w3sum,stratRaw(:,1));
tsResamp=resample(tsTest,[-0.2:0.01:1.19]);
tsDat=tsResamp.Data;
%tsDat=flip(tsDat);

tsType=timeseries(typeSum,typeStratRaw(:,1));
tsTypeResamp=resample(tsType,[-0.2:0.01:1.19]);
tsTypeDat=tsTypeResamp.Data;
%tsTypeDat=flip(tsTypeDat);
toffaList=type2cid({'rgc'},{'4ow'},curTis);

testCid=2002;
testCid=toffaList{:};
testVox=getCidVox(testCid,1,dsObj,curTis);
curVox=[[] [] []];
for y=1:length(testVox)
    curVox=vertcat(curVox,testVox{y});
end
testVox=testVox{1};
testVox=curVox;
testVox=double(testVox)./[10 10 10];
[g,n,testDepths]=getIPLdepth(testVox(:,3),testVox(:,1),testVox(:,2),[],[]);
testHistDat=histc(testDepths,[-0.2:0.01:1.19]);

dat1=testHistDat/(max(testHistDat)*1.3);
dat2=typeInterp/(max(typeInterp)*1.3);

xAxisVals=[-0.2:0.01:1.19];
figure();
hold on
%subplot(2,1,1)
plot(xAxisVals,dat1);
%xlim([0 .8])
%plot(w3sum);
%subplot(2,1,2)
plot(xAxisVals,dat2);
xlim([0 .8])

%% loop through and check some of them to make sure the fit is good.
%get the colors going correctly

plotIndividuals=1;

cellStr=struct2cell(rgcDat);
typeList=cellStr(2,1,:);
typeList=squeeze(typeList);
uniqueTypes=unique(typeList);
typeStrList={"4ow"};%,"4i","51","37","63"};
figure();
hold on
for w=1:length(typeStrList)
    curStr=typeStrList{w};
    idList=find(contains(typeList,curStr));
    typeStratRaw=[rgcDat(idList).strat];
    typeStratRaw=flip(typeStratRaw,1);
    typeStratDat=typeStratRaw(:,[2:2:end]);
    stratDat=stratRaw(:,[2:2:end]);
    if plotIndividuals
        for u=1:length(idList)
            typeSum=typeStratDat(:,u);
            typeInterp=interp1(typeStratRaw(:,1),typeSum,[-0.2:0.01:1.19]);
            curTypeDat=typeInterp/(max(typeInterp(20:100))*1.3);
            p0=plot(xAxisVals,curTypeDat);
            p0.Color=HCcolmap(w,:);
            p0.LineWidth=3;
        end
    else
        typeSum=sum(typeStratDat,2);
        typeInterp=interp1(typeStratRaw(:,1),typeSum,[-0.2:0.01:1.19]);
        curTypeDat=typeInterp/(max(typeInterp(20:100))*1.3);
        p0=plot(xAxisVals,curTypeDat);
        p0.Color=HCcolmap(w,:);
        p0.LineWidth=3;
    end
end
xlim([0 .8])
legend(typeStrList);

%% doing the same thing, but for our data
typeStrListA={"rgc","rgc","rgc","rgc","rgc"};
typeStrList={"4ow","4i","51","37","63"};
%figure();
hold on
for w=1:length(typeStrListA)
    idList=type2cid(typeStrListA{w},typeStrList{w},curTis);
    idList=idList{1};
    testVox=getCidVox(idList,1,dsObj,curTis);
    curVox=[[] [] []];
    for y=1:length(testVox)
        curVox=vertcat(curVox,testVox{y});
    end
    testVox=curVox;
    testVox=double(testVox)./[10 10 10];
    [g,n,testDepths]=getIPLdepth(testVox(:,3),testVox(:,1),testVox(:,2),[],[]);
    testHistDat=histc(testDepths,[-0.2:0.01:1.19]);
    dat1=testHistDat/(max(testHistDat(20:100))*1.3);
    p1=plot(xAxisVals,dat1,'--');
    p1.Color=HCcolmap(w,:);
    p1.LineWidth=3;
end
xlim([0 .8])
legend(typeStrList);
xline(0.47);