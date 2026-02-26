global glob





vols = glob.vol.libNames

volDir = glob.dir.Volumes;

load([volDir vols{1} '\Merge\obI.mat']);
nams = obI.nameProps.oldNames';

clear touchTag other segmentTag
for t = 1:length(nams)
    nam = nams{t};
    touchTag(t) = sum(regexp(nam,'touch')) ;
    otherTag(t) = sum(regexp(nam, 'other'));
    segmentTag(t) = sum(regexp(nam, 'Segment'));
    colorTag(t) = sum(regexp(nam, 'color'));
end


notPart = find(obI.nameProps.ofID == 0);
nams(notPart)
undefinedCids = nams(obI.nameProps.parseError>0)
rib = obI.nameProps.tag.rib;
syn = obI.nameProps.tag.syn;
junc = obI.nameProps.tag.junc;
ofID = obI.nameProps.ofID;

misc = ~(rib | syn | junc | ofID | otherTag | touchTag | segmentTag | colorTag);
badNames = nams(misc)


undefinedCids