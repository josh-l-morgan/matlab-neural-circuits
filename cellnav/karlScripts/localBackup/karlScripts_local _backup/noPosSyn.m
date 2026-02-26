allSynSegids=zeros(length(highobI.nameProps.synProp),1);
allSynPos=zeros(length(highobI.nameProps.synProp),3);
allSynNames=repmat("none",length(highobI.nameProps.synProp),1);


for i=2:length(highobI.nameProps.synProp)
    curProp=highobI.nameProps.synProp{i};
     if isfield(curProp,'segID')
        allSynSegids(i)=curProp.segID;
        allSynPos(i,:)=allColPos(curProp.segID,:);
     end
     if isfield(curProp,'name')
        allSynNames(i)=curProp.name;
     end
    
end

nopos=find(allSynPos(:,1)<1);

