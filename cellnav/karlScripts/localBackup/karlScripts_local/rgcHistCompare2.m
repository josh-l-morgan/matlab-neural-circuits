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


xAxisVals=[-0.2:0.01:1.19];

%% out of curiosity, how does the vertex distribution compare to the volume?
fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
testCid=1070;
curFv=load([fvDir num2str(testCid) '.mat']);
%curFv.fv.
% check the skeletons also
testVox=getCidVox(testCid,1,dsObj,curTis);
curVox=testVox{1};
curVox=curVox(curVox(:,1)>0&curVox(:,2)>0,:);
%testVoxSkel=horzcat(curVox(:,1).*10,curVox(:,2)*10,testDepths.*1000);
%testVoxSkel=testVoxSkel(testVoxSkel(:,1)>0&testVoxSkel(:,2)>0,:);
%figure(); scatter3(testVoxSkel(:,1),testVoxSkel(:,2),testVoxSkel(:,3));
figure(); scatter3(curVox(:,1),curVox(:,2),curVox(:,3));
emptyImg=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
emptyImg=logical(emptyImg);
indList=sub2ind([ceil(max(curVox(:,1))),ceil(max(curVox(:,2))),ceil(max(curVox(:,3)))], ...
    curVox(:,1),curVox(:,2),curVox(:,3));
indList=uint64(indList);
%for i=1:length(indList)
%    emptyImg(indList(i))=1;
%end
emptyImg(indList)=1;
emptyImg=reshape(emptyImg,ceil(max(curVox(:,1))),ceil(max(curVox(:,2))),ceil(max(curVox(:,3))));
%figure(); isosurface(emptyImg);
testSkel=bwskel(emptyImg);
testPerim=bwperim(emptyImg);
testBP=bwmorph3(testSkel,'branchpoints');


emptyImg2=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
emptyImg2(testSkel>0)=1;
skelIndList=find(testSkel>0);
%[v1,v2,v3]=ind2sub([ceil(max(testSkel(:,1))),ceil(max(testSkel(:,2))),ceil(max(testSkel(:,3)))],find(testSkel>0));
[v1,v2,v3]=ind2sub(size(testSkel),find(testSkel>0));
skelLocs=horzcat(v1,v2,v3);

emptyImg3=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
emptyImg3(testPerim>0)=1;
skelIndList=find(testPerim>0);
%[v1,v2,v3]=ind2sub([ceil(max(testSkel(:,1))),ceil(max(testSkel(:,2))),ceil(max(testSkel(:,3)))],find(testSkel>0));
[vs1,vs2,vs3]=ind2sub(size(testPerim),find(testPerim>0));
skinLocs=horzcat(vs1,vs2,vs3);

[vb1,vb2,vb3]=ind2sub(size(testBP),find(testBP>0));
bpLocs=horzcat(vb1,vb2,vb3);
%emptyImg4=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
%emptyImg4(testPerim>0)=1;
%skelIndList=find(testPerim>0);
%[v1,v2,v3]=ind2sub([ceil(max(testSkel(:,1))),ceil(max(testSkel(:,2))),ceil(max(testSkel(:,3)))],find(testSkel>0));
%[vs1,vs2,vs3]=ind2sub(size(testPerim),find(testPerim>0));
%skinLocs=horzcat(vs1,vs2,vs3);

%emptyImg3=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
%emptyImg3(testPerim>0)=1;
%skelIndList=find(testPerim>0);
%[v1,v2,v3]=ind2sub([ceil(max(testSkel(:,1))),ceil(max(testSkel(:,2))),ceil(max(testSkel(:,3)))],find(testSkel>0));
%[vs1,vs2,vs3]=ind2sub(size(testPerim),find(testPerim>0));
%skinLocs=horzcat(vs1,vs2,vs3);

[e,r,skelDepths]=getIPLdepth(skelLocs(:,3)./10,skelLocs(:,1)./10,skelLocs(:,2)./10,[],[]);
volIso=isosurface(emptyImg,0.5);
volIso2=volIso;
volIso2.vertices=volIso.vertices(:,[2 1 3]);
figure();
hold on
%t=emptyImg2(:,:,250);
%image(t);
%delTest=delaunayTriangulation(skinLocs);
%trimTest=trimesh(delTest);
%tetTest=tetramesh(delTest,'FaceAlpha',0.3);
%this sort of works, but it makes a very goofy skeleton
%volAlphShap=alphaShape(double(curVox),'HoleThreshold',5,'RegionThreshold',5);
skelScat=scatter3(v1,v2,v3,5,'c.');
bpScat=scatter3(vb1,vb2,vb3,5,'m*');
pat1=patch(volIso2);
pat1.FaceColor=[1 0 1];
pat1.LineStyle='none';
pat1.EdgeColor='none';
pat1.FaceAlpha=0.2;
%plot(volAlphShap);
%volScat=scatter3(curVox(:,1),curVox(:,2),curVox(:,3),3,'m.');
%volScat.MarkerEdgeAlpha=0.1;
%volScat.MarkerFaceAlpha=0.1;
%scatter3(v1,v2,skelDepths);
testHistDat=histc(skelDepths,[-0.2:0.01:1.19]);


%% test out some plots
if 0
typeSum=sum(typeStratDat,2);
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

figure();
hold on
%subplot(2,1,1)
plot(xAxisVals,dat1);
%xlim([0 .8])
%plot(w3sum);
%subplot(2,1,2)
plot(xAxisVals,dat2);
xlim([0 .8])
end
%% loop through and check some of them to make sure the fit is good.
%get the colors going correctly

plotIndividuals=1;

cellStr=struct2cell(rgcDat);
typeList=cellStr(2,1,:);
typeList=squeeze(typeList);
uniqueTypes=unique(typeList);
typeStrList={"37v","37d","37c","37r","5si"};
%typeStrList={"4ow","4i","51","37","63"};
%typeStrList={"4i","4ow"};
%typeStrList={""};
figure();
subplot(2,1,1)
hold on
for w=1:length(typeStrList)
    curStr=typeStrList{w};
    idList=find(contains(typeList,curStr));
    typeStratRaw=[rgcDat(idList).strat];
    typeStratRaw=flip(typeStratRaw,1);
    typeStratDat=typeStratRaw(:,[2:2:end]);
    %stratDat=stratRaw(:,[2:2:end]);
    if plotIndividuals
        for u=1:length(idList)
            typeSum=typeStratDat(:,u);
            typeInterp=interp1(typeStratRaw(:,1),typeSum,[-0.2:0.01:1.19]);
            curTypeDat=typeInterp/(max(typeInterp(20:100))*1.3);
            p0=plot(xAxisVals,curTypeDat);
            %p0.Color=HCcolmap(w,:);
            %p0.LineWidth=3;
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
xline(0.47)
legend(typeStrList);

% doing the same thing, but for our data
typeStrListA={"rgc","rgc"};%,"rgc","rgc","rgc","rgc"};
typeStrList={"5si","51"};%,"4i","51","37","63"};
%figure();
subplot(2,1,2)
hold on
method=1;
if method==1
    for w=1:length(typeStrListA)
        idList=type2cid(typeStrListA{w},typeStrList{w},curTis);
        idList=idList{1};
        testVox=getCidVox(idList,1,dsObj,curTis);
        curVox=[[] [] []];
        if plotIndividuals
            for u=1:length(idList)
                curVox=testVox{u};
                curVox=double(curVox)./[10 10 10];
                [g,n,testDepths]=getIPLdepth(curVox(:,3),curVox(:,1),curVox(:,2),[],[]);
                testHistDat=histc(testDepths,[-0.2:0.01:1.19]);
                dat1=testHistDat/(max(testHistDat(20:100))*1.3);
                p1=plot(xAxisVals,dat1,'--');
                %p1.Color=HCcolmap(w,:);
                %p1.LineWidth=3;
            end
        else
            
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
    end
elseif method==2
    for w=1:length(typeStrListA)
        idList=type2cid(typeStrListA{w},typeStrList{w},curTis);
        idList=idList{1};
        testVox=getCidVox(idList,1,dsObj,curTis);
        curVox=[[] [] []];
        if plotIndividuals
            for u=1:length(idList)
                curVox=testVox{u};
                curVox=curVox(curVox(:,1)>0&curVox(:,2)>0,:);
                emptyImg=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
                emptyImg=logical(emptyImg);
                indList=sub2ind([ceil(max(curVox(:,1))),ceil(max(curVox(:,2))),ceil(max(curVox(:,3)))], ...
                    curVox(:,1),curVox(:,2),curVox(:,3));
                indList=uint64(indList);
                %for i=1:length(indList)
                %    emptyImg(indList(i))=1;
                %end
                emptyImg(indList)=1;
                emptyImg=reshape(emptyImg,ceil(max(curVox(:,1))),ceil(max(curVox(:,2))),ceil(max(curVox(:,3))));
                %figure(); isosurface(emptyImg);
                testSkel=bwskel(emptyImg);
                emptyImg2=zeros(uint64(ceil(max(curVox(:,1))))*uint64(ceil(max(curVox(:,2))))*uint64(ceil(max(curVox(:,3)))),1);
                emptyImg2(testSkel>0)=1;
                skelIndList=find(testSkel>0);
                [v1,v2,v3]=ind2sub([ceil(max(curVox(:,1))),ceil(max(curVox(:,2))),ceil(max(curVox(:,3)))],find(testSkel>0));
                skelLocs=horzcat(v1,v2,v3);
                [e,r,skelDepths]=getIPLdepth(skelLocs(:,3)./10,skelLocs(:,1)./10,skelLocs(:,2)./10,[],[]);
                                
                testHistDat=histc(skelDepths,[-0.2:0.01:1.19]);
                dat1=testHistDat/(max(testHistDat(20:100))*1.3);
                p1=plot(xAxisVals,dat1,'--');
                %p1.Color=HCcolmap(w,:);
                %p1.LineWidth=3;
            end
        else
            
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
        
    end

end

    xlim([0 .8])
    legend(num2cell(string(idList)));
    xline(0.47);

%% Is this a function of the fit?

% going to test the 4ow
figure()
title("are Depth Estimates of 4ow branches planar?")
hold on
plotIndividuals=1;
skipFact=5;
dmap=jet(100);
boundaries={[30 45],[[23 34];[59 67]]}; 
typeStrListA={"rgc","rgc","rgc"};%,"rgc","rgc","rgc","rgc"};
typeStrList={"4ow","37","4i"};%,"4i","51","37","63"};
for w=1:length(typeStrListA)
    idList=type2cid(typeStrListA{w},typeStrList{w},curTis);
    idList=idList{1};
    testVox=getCidVox(idList,1,dsObj,curTis);
    % idList=type2cid(typeStrListA{w},typeStrList{w},curTis);
    % idList=idList{1};
    % testVox=getCidVox(idList,skipFact,dsObj,curTis);
    strat=[1 2];
    if plotIndividuals
        for u=1:length(idList)
            curVox=testVox{u};
            curVox=double(curVox)./[10 10 10];
            [g,n,testDepths]=getIPLdepth(curVox(:,3),curVox(:,1),curVox(:,2),[],[]);
            testDepths=testDepths.*100;
            %figure(); subplot(2,1,1); histogram(botbottomCols); subplot(2,1,2); histogram(topbottomCols);
            %figure(); histogram(testDepths);
            if strat(w)==1
                curBounds=boundaries{w};
                bottomCols=-curBounds(1,1)+testDepths;
                bottomCols(bottomCols>0)=0;
                topCols=testDepths-curBounds(1,2);
                topCols(topCols<0)=0;
                colVals=bottomCols+topCols;
                colVals=(colVals.*10)+50;
                colVals(colVals<1)=1;
                colVals(colVals>100)=100;
                colVals=round(colVals);
            elseif strat(w)==2
                curBounds=boundaries{w};
                botbottomCols=-curBounds(1,1)+testDepths;
                botbottomCols(botbottomCols>0)=0;
                topbottomCols=-curBounds(2,1)+testDepths;
                topbottomCols(topbottomCols>0)=0;
                topbottomCols(topbottomCols<(-curBounds(2,1)+47))=0;
                
                bottopCols=-curBounds(1,2)+testDepths;
                bottopCols(bottopCols<0)=0;
                bottopCols(bottopCols>(-curBounds(1,2)+47))=0;
                toptopCols=-curBounds(2,2)+testDepths;
                toptopCols(toptopCols<0)=0;
                
%                 figure();
%                 hold on
%                 subplot(4,1,1);
%                 histogram(botbottomCols);
%                 subplot(4,1,2);
%                 histogram(topbottomCols);
%                 subplot(4,1,3);
%                 histogram(bottopCols);
%                 subplot(4,1,4);
%                 histogram(toptopCols);
                
                colValDat=horzcat(botbottomCols,topbottomCols,bottopCols,toptopCols);
                colVals=sum(colValDat,2);
                colVals=(colVals.*3)+50;
                colVals(colVals<1)=1;
                colVals(colVals>100)=100;
                colVals=round(colVals);
                
                
            end
%             colDepths=round((testDepths).*1000);
%             colDepths(colDepths<201)=201;
%             colDepths(colDepths>600)=600;
%             colDepths=colDepths-200;
            colMat=dmap(colVals,:);
            scatter3(curVox(:,1),curVox(:,2),testDepths,3,colMat,'filled');
        end
    end
end
legend(string(idList))
xlim([70 210])
ylim([50 170])
view(90,90);
    