%functest
bpcLists=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc'},{'bc3a','bc3b','bc4','bc5o','bc5i','bc5t'},curTis);

rgcLists=type2cid({'rgc','rgc','rgc','rgc','rgc'},{'4ow','5ti','37','6sw','63'},curTis);

cidCellArrayA=cell(2,1);
cidCellArrayA{1}=vertcat(bpcLists{1},bpcLists{2},bpcLists{3});
cidCellArrayA{2}=vertcat(bpcLists{4},bpcLists{5},bpcLists{6});
axisNamesA={'OFF bpc','ON bpc'};

cidCellArrayB=cell(3,1);
cidCellArrayB{1}=rgcLists{1};
cidCellArrayB{2}=rgcLists{2};
cidCellArrayB{3}=rgcLists{3};
cidCellArrayB{4}=rgcLists{4};
cidCellArrayB{5}=rgcLists{5};
axisNamesB={'4ow','5ti','37','6sw','63'};

useSkels=[1 2 3 5 6];
for i=1:length(useSkels)
    curSkel=useSkels(i);
    SMcellArray{i}=allSkels{curSkel}.sm;
end

[outputDat stats]=getDRG(cidCellArrayA,cidCellArrayB,SMcellArray,[0:1:30]);

sumDat=sum(outputDat,4);

size(sumDat)

f1=figure();
hold on
%tl1=tiledlayout(size(sumDat,2),size(sumDat,3));

ylims={[0 1], [0 10]};
for i=1:size(sumDat,2)
    for j=1:size(sumDat,3)
    subplot(1,size(sumDat,3),j)
    plot((sumDat(:,i,j))/stats(j)*2)
    hold on
    %title([axisNamesA{i} ' x ' axisNamesB{j}])
    title(['BPC polarity x ' axisNamesB{j}])
    ylim(ylims{1});
    end
end
    



%% Next section
%what are the BPC inputs of the different RGC types?
%important cells are 2002, 3051, and 3119 (4ow,5ti,37)

%find the bpcs preVG
allVGcids=[2 3 4 5 10 11 13 14 20];
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
b2vIDs=find(allPreTypes{1}'==7&ismember(curTis.syn.edges(:,1),allVGcids));
b2vCids=curTis.syn.edges(b2vIDs,2);
b2vCids=unique(b2vCids);

cidList=[2002 3051 3119];
for k=1:length(cidList)
    curCid=cidList(k);
    curInputs=find(allPreTypes{1}'==7&curTis.syn.edges(:,1)==curCid);
    preBPCcids=curTis.syn.edges(curInputs,2);
    knownBPCbool=ismember(preBPCcids,b2vCids);
    mean(knownBPCbool)
end