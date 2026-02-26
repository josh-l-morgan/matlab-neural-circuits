function outTis=findBadSynapse(curTis)
%global tis
%curTis=tis;
keepGoing=1;

while keepGoing
    %% find synapses with no position
    badSyn=find(curTis.syn.pos(:,1)<1);
    badSynPos=curTis.syn.pos(badSyn,:);
    badSynObID=curTis.syn.obID(badSyn);
    badSynName=curTis.obI.nameProps.oldNames(badSynObID);
    badSynSeg=obI.fuse.obSource(badSynObID);
    badSynSegName=obI.fuse.exportDir(badSynSeg);
    
    %% find synapses with identical positions
    allSynPos=curTis.syn.pos;
    distMat=pdist2(allSynPos,allSynPos);
    noDistSyns=find(distMat<0.2);
    self=sub2ind(size(distMat),[1:size(distMat,1)],[1:size(distMat,2)]);
    noDistSynsClean=noDistSyns(~ismember(noDistSyns,self));
    [badx,bady]=ind2sub(size(distMat),noDistSynsClean);
    badInds=horzcat(badx,bady);
    %% debugging
    debug=1;
    if debug
        absIndDif=abs(badx-bady);
        checkInds=badInds(find(absIndDif>1),:);
        for checkIt=1:size(checkInds,1)
            checkPair=checkInds(checkIt,:);
            
            
            
        end
    end
    %% combine the two lists
    realBadInds=[];
    badIndsCln=badInds(~ismember(badInds(:,2),badSyn),:);
    i=1;
    while i<=size(badIndsCln,1)
        curCompSynIDs=badIndsCln(i,:);
        badEdges=curTis.syn.edges(curCompSynIDs,:);
        badEdgesAC=aliasChecker(badEdges(:,[1 2]),aliasList);
        badPos=curTis.syn.pos(curCompSynIDs,:);
        badPosVast=badPos.*[250 250 25];
        badNames=curTis.obI.nameProps.names(curTis.syn.obID(curCompSynIDs));
        badSegs=obI.fuse.obSource(badEdges(:,3));
        badSegNames=obI.fuse.exportDir(badSegs);
        output=horzcat(badNames(:),badSegNames(:))';
        output2=[badEdgesAC(1,:) 1 badEdgesAC(2,:) 1 ;...
            badPosVast(1,:) badPosVast(2,:)];
        
        output
        output2
        pause
        if sum(badEdges(1,1:2)==badEdges(2,1:2))==2
            realBadInds=[realBadInds;curCompSynIDs(1)];%;curCompSynIDs(2)];
            i=i+1;
            %a=curTis.syn.pos(curComSynIDs(1),[2 1 3])
            %l=curTis.syn.edges(curCompSynIDs(1),:)
            %a=curTis.syn.pos(curCompSynIDs(1),[2 1 3])
            %b=a.*[250 250 25];
            %clipboard('copy',b);
            %pause;
        end
        i=i+1;
    end
    
    %% Need to update to look for ribbons that have the same
    % same bpc
    % close location
    % at least 1 same target
    % remove the one with less / combine the downstream cells
    
    %% update the current tis file and re-run
    
    noPosSyns=badSyn;
    realBadInds=[realBadInds;noPosSyns];
    removeInds=unique(realBadInds);
    
    if size(removeInds,1)==0
        keepGoing=0;
    end
    keepInds=[1:length(curTis.syn.pos(:,1))];
    keepInds=keepInds(~ismember(keepInds,removeInds));
    curTis.syn.edges=curTis.syn.edges(keepInds,:);
    curTis.syn.pre=curTis.syn.pre(keepInds,:);
    curTis.syn.post=curTis.syn.pos(keepInds,:);
    curTis.syn.obID=curTis.syn.obID(keepInds,:);
    curTis.syn.synType=curTis.syn.synType(keepInds,:);
    curTis.syn.preClass=curTis.syn.preClass(keepInds,:);
    curTis.syn.postClass=curTis.syn.postClass(keepInds,:);
    curTis.syn.pos=curTis.syn.pos(keepInds,:);
    curTis.syn.synPosRaw=curTis.syn.synPosRaw(keepInds,:);
    curTis.syn.synPosDS=curTis.syn.synPosDS(keepInds,:);
        
end

outTis=curTis;
end