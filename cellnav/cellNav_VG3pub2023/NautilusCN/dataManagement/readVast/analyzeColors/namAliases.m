function[newNam,newCids,cids] = namAliases(nam,aliases)

%%Find cids and replace cids with proper aliases creating new nam and new
%%cid list

if ~exist('cidTag','var')
    cidTag = 'cid';
end
if ~exist('allowSpace','var')
    allowSpace = 1;
end

%% Get all cids in nam
cidPos = regexp(nam,cidTag);
cids = zeros(length(cidPos),1);
parseError = 0;
isSpace = zeros(size(cids));
cidStart = cidPos * 0;
for r = 1:length(cidPos)

    postCid = nam(cidPos(r)+3:end);
    cidN = sscanf(postCid,'%d');
    if ~isempty(cidN)
        cids(r) = cidN(1);
        starts = regexp(postCid,num2str(cidN(1)));
        cidStart(r) = cidPos(r)+2 + starts(1);
    end
    
    if ~isempty(cidN)
        cids(r) = cidN(1);
    else
        cids(r) = 0;
        parseError = 1;
    end
end

%% find replacement cids
newCids = cids;
for i = 1:length(aliases)
    alias = aliases{i};
    for a = 2:length(alias)
        newCids(newCids==alias(a)) = alias(1);
    end
end


%% make modified nam
newNam = nam;
L1 = length(nam);
for i = 1:length(newCids)
    L2 =length(newNam);
    Ld = L2 - L1;
    beforeS = newNam(1:cidStart(i)-1 + Ld);
    cidS = num2str(cids(i));
    afterS = newNam(cidStart(i)+length(cidS)+Ld:end);
    newNam = [beforeS num2str(newCids(i)) afterS];
end










