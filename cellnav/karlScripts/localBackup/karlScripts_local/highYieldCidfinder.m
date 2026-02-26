%% Look by connectivity
vgcAllCids=[2 3 4 5 10 11 13 14 20];
allTargCids=find(ismember(curTis.syn.edges(:,2),vgcAllCids));
allTargType=cid2type(curTis.syn.edges(:,1),curTis);
allCids=find(allTargType{1}==1 | allTargType{1}==0 | allTargType{1}==8);
allv2r=intersect(allTargCids,allCids);
allTargRGCCids=curTis.syn.edges(allv2r,1);
test=unique(allTargRGCCids);
targRGClist=test;
targRGCcounts=targRGClist;
for i=1:length(targRGClist)
    curRGCcid=targRGClist(i);
    curRGCcount=sum(allTargRGCCids==curRGCcid);
    targRGCcounts(i)=curRGCcount;
end
outputType=cid2type(targRGClist,curTis);
outputDat=horzcat(targRGClist,targRGCcounts,outputType{1}',outputType{3}');
outputDatSrtd=sortrows(outputDat,[2 3]);

types=zeros(length(unique(outputDat(:,3))),2);
types(:,1)=unique(outputDat(:,3));
for k=1:size(types,1)
    types(k,2)=sum(outputDat(outputDat(:,3)==types(k,1),2));
end





%% Look by Size
load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat');
allRGCcids=type2cid({'rgc'},{'all'},curTis);
%allUNKcids=type2cid({'unk',' '},{'all','all'},curTis);
allRGCcids=allRGCcids{1};
allRGCcids=outputDatSrtd(:,1);
%allOtherCids=allUNKcids{:};
allRGCvox=getCidVox(allRGCcids,1,dsObj,curTis);
voxCounts=zeros(size(allRGCcids));
for i=1:length(allRGCvox)
    voxCounts(i)=length(allRGCvox{i});
end
cidsCounts=horzcat(allRGCcids,voxCounts);
[s,d]=sort(voxCounts,'descend');
cidCntSrt=cidsCounts(d,:);
subtypes=cid2type(cidCntSrt(:,1),curTis);
subText=subtypes{4}';
subtypes=subtypes{3}';

cidCntTyp=horzcat(cidCntSrt,subtypes);


%% make a plot
ODf=horzcat(outputDatSrtd,voxCounts);
dx=0.05;
dy=0.05;
cmp=turbo(16);
f11=figure();
hold on;
scatter(ODf(:,2),ODf(:,5),10,cmp(ODf(:,3)+1));
text(ODf(:,2)+dx, ODf(:,5)+dy, num2str(ODf(:,1)), 'Fontsize', 10);






%%
while 0
    %skip next section normally
%%

dsObj=highObj;
curTis=highTis;
bins=[0:15:300];

%finding the others
allTypeDat=cid2type(curTis.cells.cids,curTis);
rgcCids=curTis.cells.cids(find(allTypeDat{1}==1));
otherCids=curTis.cells.cids(find(allTypeDat{1}==0 | allTypeDat{1}==8));
otherCids=otherCids(~ismember(otherCids,vgcAllCids));

allCids=horzcat(rgcCids,otherCids);
%get the sizes of all the cells
allVox=getCidVox(allCids,1,dsObj,curTis);
allHist=zeros(length(allVox),length(bins)-1);
voxCounts=zeros(size(allCids));
for i=1:length(allVox)
    voxCounts(i)=length(allVox{i});
    voxDat=allVox{i};
    curHist=histcounts(voxDat(:,3),bins);
    allHist(i,:)=curHist;
end
cidsCounts=vertcat(allCids,voxCounts);
cidsCounts=cidsCounts';
[s,d]=sort(voxCounts,'descend');
cidCntSrt=cidsCounts(d,:);
VinputNum=getVGinputs(cidCntSrt(:,1),curTis);
ttype=cid2type(cidCntSrt(:,1),curTis);
typeVal=ttype{1};
%largerThanTenK=cidCntSrt(cidCntSrt(:,2)>10000,1);

largDat=horzcat(cidCntSrt,VinputNum,typeVal');

figure();
%hold on;
for i=1:length(allCids)
    plot(bins(2:end)-1,allHist(i,:));
    xticks(bins(2:end)-1);
    xticklabels(bins(2:end));
    title(num2str(allCids(i)));
    pause();
end

%make a plot
dx=0.05;
dy=0.05;
cmp=turbo(16);
f10=figure();
hold on;
scatter(largDat(:,2),largDat(:,3),10,cmp(largDat(:,4)+1));
text(largDat(:,2)+dx, largDat(:,3)+dy, num2str(largDat(:,1)), 'Fontsize', 10);


%find everyone with a traced soma in the INL

allVox=getCidVox(allCids,1,dsObj,curTis);
allHistDat=1;





%get a plane for the overall volume
testpts=zeros(1,3);
morePoints={};
morePoints{1}=testpts;
morePoints{2}=testpts;
morePoints{3}=testpts;
morePoints{4}=testpts;
morePoints{5}=testpts;

results={};
for i=1:length(morePoints)
testpts=morePoints{i};
testptsTd=testpts(:,[3 2 1])./[25 250 250];
x=testptsTd(:,1);
y=testptsTd(:,2);
z=testptsTd(:,3);
DM = [x, y, ones(size(z))];                             % Design Matrix
B = DM\z;                                               % Estimate Parameters
[X,Y] = meshgrid(linspace(min(x),max(x),50), linspace(min(y),max(y),50));
Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));
results{i,1}=B;
results{i,2}=Z;
end

curPlane=[-1 45.0476388697178 0.969913040048901 -1121.00155821979];
[X,Y] = meshgrid(linspace(0,200,10), linspace(0,200,10));
zPlane=(-double(x)*curPlane(2)-double(y)*curPlane(3) ...
    -curPlane(4))/curPlane(1);

%test the planes
pltpts=morePoints{2};
pltpts=pltpts(:,[3 2 1])./[25 250 250];
f11=figure();
hold on;
compareMorph(f11,[1110],[]);
scatter3(pltpts(:,1),pltpts(:,2),pltpts(:,3),50,'ro');

%%
end
allTypeDat=cid2type(curTis.cells.cids,curTis);
allRGCcids2=curTis.cells.cids(find(allTypeDat{1}==1));