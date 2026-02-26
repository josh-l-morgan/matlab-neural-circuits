
HCcolmap=load('HCcolmap.mat');
    HCcolmap=HCcolmap.HCcolmap;
    HCcolmap=HCcolmap./255;

%% get the list of the synapses
bpcLists=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc'},{'bc3a','bc3b','bc4','bc5o','bc5i','bc5t','bc6','bcoff','bcon'},curTis);

rgcLists=type2cid({'rgc','rgc','rgc','rgc','rgc','rgc'},{'4ow','4i','5ti','37','6sw','63'},curTis);



cidCellArrayA={};
% cidCellArrayA{1}=vertcat(bpcLists{1},bpcLists{2},bpcLists{3},bpcLists{8});
% cidCellArrayA{2}=vertcat(bpcLists{4},bpcLists{5},bpcLists{6},bpcLists{7},bpcLists{9});
% cidCellArrayA{3}=0;
% axisNamesA={'OFF bpc','ON bpc','amc in'};

cidCellArrayA=bpcLists;
axisNamesA={'bc3a','bc3b','bc4','bc5o','bc5i','bc5t','bc6','bcoff','bcon'};


cidCellArrayB={};
% cidCellArrayB{1}=rgcLists{1};
% cidCellArrayB{2}=rgcLists{2};
% cidCellArrayB{3}=rgcLists{3};
% cidCellArrayB{4}=rgcLists{4};
% cidCellArrayB{5}=rgcLists{5};
% cidCellArrayB{6}=rgcLists{6};
cidCellArrayB=rgcLists;
cidCellArrayB{length(cidCellArrayB)+1}=0;
axisNamesB={'4ow','4i','5ti','37','6sw','63','amc'};

useSkels=[1 2 3 5 6];
for i=1:length(useSkels)
    curSkel=useSkels(i);
    SMcellArray{i}=allSkels{curSkel}.sm;
end

binEdges=[0:.5:29.5];

[outputDat f1]=getDRG2(cidCellArrayA,axisNamesA,cidCellArrayB,axisNamesB,SMcellArray,binEdges,3);

%% parse that output
totalG2G=zeros(size(outputDat(1).g2g{1},1),size(outputDat(1).g2g{1},2),length(cidCellArrayB));
totalCounts=zeros(length(cidCellArrayB),1);
for i=1:length(useSkels)
   for j=1:length(cidCellArrayB)
   totalG2G(:,:,j)=totalG2G(:,:,j)+outputDat(i).g2g{j};
   totalCounts(j)=totalCounts(j)+outputDat(i).counts{j};
   end
end

fab=figure();
hold on
sgtitle('Neighborhood BPC input density for outputs to RGC subtypes')
for i=1:length(cidCellArrayB)
subplot(2,4,i)
%subplot(6,6,30+i)
hold on
for j=1:length(cidCellArrayA)%[1 2 3 5]
    plotz=plot(binEdges,totalG2G(j,:,i)./totalCounts(i),'LineWidth',2);
    %yline(means(i)/100);
    if ~ismember(j,[1,2,3,8])
        plotz.LineStyle='--';
    end
end
title([axisNamesB{i} ' n=' num2str(totalCounts(i))]);
% for i=1:5%[1 2 3 5]
%     %plot(binEdges,test(i,:)./100,'LineWidth',2);
%     yline(means(i)/100,':');
% end
xlabel('distance (um)');
ylabel('density (#syn/um)');
xlim([0 10])
ylim([0 .008])
if i==length(cidCellArrayB)
    %subplot(2,4,8);
    leggy=legend(axisNamesA)
end
end
leggy.Position(1)=0.8;
%legend({'bpcIn','amcIn','rgcOut','amcOut'})

%% get the histogram overlap
%getHisto
%getHist
histBins=[-.20:.01:1.20];
histDatA=zeros(length(histBins)-1,length(cidCellArrayA));
histDatB=zeros(length(histBins)-1,length(cidCellArrayB));
for m=1:length(cidCellArrayA)
    curCidList=cidCellArrayA{m};
    listVox=getCidVox(curCidList,1,dsObj,curTis);
    curTypeHistDat=zeros(length(histBins)-1,1);
    for n=1:length(listVox)
        curVox=listVox{n};
        curVox=double(curVox)/10;
        [g l curCidDepths]=getIPLdepth(curVox(:,3),curVox(:,1),curVox(:,2),[],[]);
        curCidHistDat=histcounts(curCidDepths,histBins)';
        curTypeHistDat=curTypeHistDat+curCidHistDat;
    end
    histDatA(:,m)=curTypeHistDat;
end

for m=1:length(cidCellArrayB)
    curCidList=cidCellArrayB{m};
    listVox=getCidVox(curCidList,1,dsObj,curTis);
    curTypeHistDat=zeros(length(histBins)-1,1);
    for n=1:length(listVox)
        curVox=listVox{n};
        curVox=double(curVox)/10;
        [g l curCidDepths]=getIPLdepth(curVox(:,3),curVox(:,1),curVox(:,2),[],[]);
        curCidHistDat=histcounts(curCidDepths,histBins)';
        curTypeHistDat=curTypeHistDat+curCidHistDat;
    end
    histDatB(:,m)=curTypeHistDat;
end

testFig=figure();
hold on
for i=1:length(cidCellArrayA)
    p1=plot(histBins(1:end-1),histDatA(:,i));
    p1.LineWidth=2;
end

for i=1:length(cidCellArrayB)
    p2=plot(histBins(1:end-1),-histDatB(:,i),'--');
    p2.LineWidth=2;
end
legend({axisNamesA{:} axisNamesB{:}},'Location','east');
xlim([0.15 0.85]);

    
%% see what the eyewire arbor overlaps look like.
eyewireNames=axisNamesB(1:end-1);
for k=1:length(eyewireNames)
    curEyeType=eyewireNames{k};
end

%% go through the outputs and plot the input neighborhoods
while 0
plotz=1;
synRad=10;
for i=1:length(SMcellArray)
    curSM=SMcellArray{i};
    curCidFv=curSM.nep.fv;
    for m=1:size(outputDat(1).locsB,2)
        outputIDs=outputDat(i).idsB{m};
        f1=figure();
        for j=1:length(outputIDs)
            hold off
            hold on
            curTipID=outputIDs(j);
            tipPos=outputDat(i).locsB{m}(j,:);
            dists2syns=curSM.syn2Skel.syn2SynDist(:,curTipID);
            %dists2synsA=dists2syns(inputIDs);
            %closeAids=find(dists2synsA<synRad);
            %figure(); histogram(dists2syns,[0:150]);
            %tipPos=curSM.nep.pos(curTipID,:);
            %tipZ=nepDepths(curTipID,:);
            figure(f1);
            scatter3(tipPos(1),tipPos(2),tipPos(3),250,'cp','filled');
            plotRad=10;
            xlim([tipPos(1)-plotRad tipPos(1)+plotRad ])
            ylim([tipPos(2)-plotRad tipPos(2)+plotRad ])
            zlim([tipPos(3)-plotRad tipPos(3)+plotRad ])
            axis vis3d
            fvPatch=patch(curCidFv);
            fvPatch.EdgeColor='none';
            fvPatch.FaceColor=[.2 .2 .2];
            fvPatch.FaceAlpha=.1;
            %find all nearby synapses
            %closeSyns=find(dists2syns<synRad);
            %get positions of close syns
            % closeSynPos=curSM.syn.pos(closeSyns,:);
            %scatter3(inputLocs(closeAids,1),inputLocs(closeAids,2),inputLocs(closeAids,3),100,'mo','filled');
            acols=[[1 0 1];[1 .5 0]];
            for n=1:size(outputDat(1).locsA,2)
                inputIDs=outputDat(i).idsA{n};
                inputLocs=outputDat(i).locsA{n};
                dists2synsA=dists2syns(inputIDs);
                closeAids=find(dists2synsA<synRad);
                scatter3(inputLocs(closeAids,1),inputLocs(closeAids,2),inputLocs(closeAids,3),100,acols(n,:),'filled');
            end
            while 0
                if 0
                    %find the kinds of synapses so we can plot them correctly
                    closeIn=intersect(curSM.syn.isIn,closeSyns);
                    closeOut=intersect(curSM.syn.isOut,closeSyns);
                    closeOutTypes=postTypes{1}(closeOut);
                    closeInTypes=preTypes{1}(closeIn);
                    
                    
                    amcIn=closeIn(closeInTypes==8|closeInTypes==0);
                    bpcIn=closeIn(closeInTypes==7);
                    amcOut=closeOut(closeOutTypes==8|closeOutTypes==0);
                    rgcOut=closeOut(closeOutTypes==1);
                    if plotz2
                        figure(f2);
                        if ~isempty(amcIn)
                            scatter(dists2syns(amcIn),repmat(j,[length(dists2syns(amcIn)) 1]), ...
                                20,'mv','filled');
                        end
                        if ~isempty(bpcIn)
                            scatter(dists2syns(bpcIn),repmat(j,[length(dists2syns(bpcIn)) 1]), ...
                                20,'cv','filled');
                        end
                        if ~isempty(amcOut)
                            scatter(dists2syns(amcOut),repmat(j,[length(dists2syns(amcOut)) 1]), ...
                                20,'m^','filled');
                        end
                        if ~isempty(rgcOut)
                            scatter(dists2syns(rgcOut),repmat(j,[length(dists2syns(rgcOut)) 1]), ...
                                20,'g^','filled');
                        end
                        
                    end
                    if plotz
                        
                        synScale=75;
                        scaleSyns=0;
                        synSizes=[synScale synScale synScale synScale];
                        if scaleSyns
                            synSizes(1)=synScale./dists2syns(amcIn);
                            synSizes(2)=synScale./dists2syns(bpcIn);
                            synSizes(3)=synScale./dists2syns(amcOut);
                            synSizes(4)=synScale./dists2syns(rgcOut);
                        end
                        
                        closeSynCols=repmat([0 0 0],[length(closeSyns) 1]);
                        closeSynCols(ismember(closeSyns,closeIn),:)=repmat(HCcolmap(4,:),[length(closeIn) 1]);
                        closeSynCols(ismember(closeSyns,closeOut),:)=repmat(HCcolmap(9,:),[length(closeOut) 1]);
                        %plot'em
                        %closeSynScat=scatter3(closeSynPos(:,1),closeSynPos(:,2),closeSynPos(:,3),100./dists2syns(closeSyns),closeSynCols,'o','filled');
                        closeSynScat=scatter3(curSM.syn.pos(amcIn,1),curSM.syn.pos(amcIn,2),curSM.syn.pos(amcIn,3),synSizes(1),'mv','filled');
                        closeSynScat=scatter3(curSM.syn.pos(bpcIn,1),curSM.syn.pos(bpcIn,2),curSM.syn.pos(bpcIn,3),synSizes(2),'cv','filled');
                        closeSynScat=scatter3(curSM.syn.pos(amcOut,1),curSM.syn.pos(amcOut,2),curSM.syn.pos(amcOut,3),synSizes(3),'m^','filled');
                        closeSynScat=scatter3(curSM.syn.pos(rgcOut,1),curSM.syn.pos(rgcOut,2),curSM.syn.pos(rgcOut,3),synSizes(4),'g^','filled');
                        %closeSynScat.Marker=cell(size(closeSynPos,1),1);
                        %           closeSynScat.CData=
                        %             for l=1:size(closeSynPos,1)
                        %                 curSynPos=closeSynPos(l,:);
                        %                 p=plot3([curSynPos(3) tipPos(3)],[curSynPos(1) tipPos(1)], ...
                        %                     [curSynPos(2) tipPos(2)],'--','LineWidth',25/dists2syns(closeSyns(l)));
                        %                 p.LineWidth
                        %             end
                    end
                end
            end
            pause();
        end
    end
    %title(num2str(curCid))
end
end