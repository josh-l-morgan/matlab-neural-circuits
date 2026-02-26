function[cids cidPos parseError] = getCids(nam,cidTag,allowSpace)

if ~exist('cidTag','var')
    cidTag = 'cid';
end
if ~exist('allowSpace','var')
    allowSpace = 1;
end


cidPos = regexp(nam,cidTag);
cids = zeros(length(cidPos),1);
parseError = 0;

for r = 1:length(cidPos)
    
    cidN = sscanf(nam(cidPos(r)+3:end),'%d');
    if ~isempty(cidN)
        cids(r) = cidN(1);
    elseif allowSpace
        cidN = sscanf(nam(cidPos(r)+4:end),'%d');
    end
    if ~isempty(cidN)
        cids(r) = cidN(1);
    else
        cids(r) = 0;
        parseError = 1;
    end
end