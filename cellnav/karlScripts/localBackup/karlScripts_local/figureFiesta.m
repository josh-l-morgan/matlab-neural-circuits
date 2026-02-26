%% A pretty graph of bipolar input synapse depth distributions
% parameters
depthStep=0.01;
smoothBool=1;

allVG3cids=[2 3 4 5 6 10 11 13 14 20];
allSynPos=curTis.syn.pos;
[n,m,allSynDepth]=getIPLdepth(allSynPos(:,3),allSynPos(:,1),allSynPos(:,2),[],[]);
bpcTypeSets={[[7 7 7 7];[ 3 4 5 19]],[[7 7 7 7 7 7];[6 7 8 19 20 21]]};
bpcOnSubtypes=[6:12,14,17,21:22,26];
bpcOffSubtypes=[1:5,15,18:20,23];
bpcOnInInds=getClassConn(7,bpcOnSubtypes,8,1,curTis);
bpcOffInInds=getClassConn(7,bpcOffSubtypes,8,1,curTis);

depthHistBins=[0:depthStep:1];

bpcOnInDepths=allSynDepth(bpcOnInInds);
bpcOffInDepths=allSynDepth(bpcOffInInds);

bpcOnHistDat=histcounts(bpcOnInDepths,depthHistBins);
bpcOffHistDat=histcounts(bpcOffInDepths,depthHistBins);

groupHistDat={bpcOnHistDat,bpcOffHistDat};

% make the actual figure
f01=figure();
hold on
bars={};
for i=1:length(groupHistDat)
    curHistDat=groupHistDat{i};
    if smoothBool
        curHistDat=smooth(curHistDat);
    end
    xaxvals=[depthStep/2:depthStep:1-depthStep/2];
    curB=bar(xaxvals*100,curHistDat);
    curB.FaceAlpha=0.4;
    bars{i}=curB;
end
%xline(47);
legend({'ON bipolar inputs','OFF bipolar inputs'});
ylabel('# of input synapses');
xlabel('IPL depth %');
title('IPL depth of bipolar inputs to VG3 arbor by polarity');
xlim([20 70])


%% old code
if 0
%old stuff
%indices
%super sketch, but I'm just saying that input&~bpc = amc
bpcSynIDs=find(preTypeDat{1}==7 & inputs'==1);
amcSynIDs=find(inputs'==1 & preTypeDat{1}~=7);
onBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcONsubs));
offBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcOFFsubs));
bpcIDs={};
bpcHistDat={};
for i=1:length(bpcSubs)
    curSub=bpcSubs(i);
    bpcIDs{i}=find(preTypeDat{1}==7 & inputs'==1 & preTypeDat{3}==curSub);
    curHist=histcounts(allSynDepth(bpcIDs{i}),depthHistBins);
    bpcHistDat{i}=curHist;
end
rgcIDs={};
rgcHistDat={};
for i=1:length(rgcSubs)
    curSub=rgcSubs(i);
    rgcIDs{i}=find(postTypeDat{1}==1 & outputs'==1 & postTypeDat{3}==curSub);
    curHist=histcounts(allSynDepth(rgcIDs{i}),depthHistBins);
    rgcHistDat{i}=curHist;
end

%%
% try to figure out the 'distances' between input and output synapse depth
% distributions
histDistMat=zeros(length(bpcHistDat),length(rgcHistDat));
for j=1:length(bpcHistDat)
    for k=1:length(rgcHistDat)
        histDistMat(j,k)=pdist2(bpcHistDat{j},rgcHistDat{k});
    end
end

%% depth dist plot
f01=figure();
hold on
subplot(2,1,1);
hold on
title('BPC input IPL depths');
for i=1:length(bpcIDs)
    curHistDat=bpcHistDat{i};
    a=area(curHistDat,'LineStyle','--');
    a.FaceAlpha = 0.2;
    %plot(1:length(curHistDat),curHistDat,'--');
end
xticks([0.5:2:19.5]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
set(gca, 'YDir','reverse');
legend(bpcSubNames);
xlim([4,15]);
%xline(9.5);
xline(10,'HandleVisibility','off');
xlabel('IPL depth');
ylabel('synapse #');

subplot(2,1,2);
hold on
title('RGC output IPL depths');
for i=1:length(rgcIDs)
    curHistDat=rgcHistDat{i};
    %plot(1:length(curHistDat),curHistDat,':');
    a=area(curHistDat,'LineStyle','--');
    a.FaceAlpha = 0.2;
end
legend(rgcSubNames);
xticks([0.5:2:19.5]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
xline(10,'HandleVisibility','off');
xlim([4,15]);
xlabel('IPL depth');
ylabel('synapse #');
%legend({bpcSubNames{:},rgcSubNames{:}});
end