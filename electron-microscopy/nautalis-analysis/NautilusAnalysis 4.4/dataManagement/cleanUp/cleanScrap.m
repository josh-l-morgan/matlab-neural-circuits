
    load([MPN 'obI.mat'])


names = cat(1,obI.colStruc.names(:))
targ =  '486'
isTarg = [];
c = 0;
for i = 1:length(names)
   nam = names{i};
    if sum(regexp(nam,targ))
        c = c + 1;
        isTarg = [isTarg i];
        foundNames{c} = nam;
    end
    
end
    


foundNames
source = obI.fuse.obSource(isTarg);
obI.fuse.files(source).subTag


%%

unique(obI.nameProps.cellNum)'

%% allIDs
sort(unique([obI.nameProps.edges(:); obI.cell.name(:)]))


%% check for no ID 

synIDs = obI.nameProps.edges;

synNum = size(synIDs,1);
IDs = obI.cell.name;

c = 0;
clear missName missSource
for i = 1:synNum
    hit = sum(IDs == synIDs(i,2));
    if ~hit
        c = c+1;
       missName(c) = synIDs(i,2);
       obId = synIDs(i,3);
       missSource{c} = obI.fuse.files(obI.fuse.obSource(obId)).subTag;
    end
    
end


% find sources 

targ = 208;

isTarg = find(obI.nameProps.cellNum == targ);
names = obI.nameProps.names(isTarg);

source = obI.fuse.obSource(isTarg);
obI.fuse.files(source).subTag
obI.fuse.files.subTag



%% get indicies

names = obI.colStruc.names;
tag = '393 part 125';
c = 0;
clear isName isSeg anch
load([MPN 'dsObj.mat'])
for i = 1:length(names)
    nam = names{i};
    if sum(regexp(nam,tag));
     c = c+1;
     isName{c} = nam;
     isSeg(c) = i;
     anch(c,:) = obI.colStruc.anchors(i,:);
        
    end
    
end

sub = dsObj(isSeg(2)).subs;

%% check for synapse

pre = 125;
post = 494;

sum((allEdges(:,1)==post) & (allEdges(:,2) == pre))
obI.fuse.files.subTag


targ = 'cell 492';

nams = obI.colStruc.names;
for i = 1:length(nams)
    nam = nams{i};
    if sum(regexp(nam,targ))
        nam
    end
end



%% tcrs missing from segmentation


synPos = getSynPos(1);
mot = getMotifs(obI);
synStruct = getSynMat(obI);
%anaMot = analyzeMotifs(obI);
taggedTCRs = synPos.tcrs;


from125 = synStruct.pre == 125;
tracedTCRs = synStruct.post(from125);

setdiff(taggedTCRs, tracedTCRs)

setdiff(tracedTCRs, taggedTCRs)
























