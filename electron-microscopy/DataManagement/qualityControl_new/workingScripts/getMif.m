function[mif] = getMif(TPN)

if ~exist('TPN','var')
    TPN = GetMyDir;
end

APN = findFolders(TPN);
allNams = {}; allTimes = []; allFolders = {};
for f = 1: length(APN)
    sprintf('finding tiles in folder %d of %d',f,length(APN))
    [aNams aTimes]= getTiles(APN{f});
    fID = ones(length(aNams),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
    allNams = cat(1,allNams, aNams);
    allTimes = cat(1,allTimes, aTimes);
end


mif.dir = TPN;

for i = 1:length(allNams)
    %pause(.1)
    nam = allNams{i};
    rs = regexp(nam,'_r');
    us= regexp(nam,'_');
    cs = regexp(nam,'-c');
    ws = regexp(nam,'_w');
    ss = regexp(nam,'_sec');
    ts = regexp(nam,'.tif');
    if ~isempty(rs) & ~isempty(cs)
        r(i) = str2num(nam(rs(1)+2:cs(1)-1));
        c(i) = str2num(nam(cs(1)+2:us(2)-1));
        w(i) = str2num(nam(ws(1)+2:ss(1)-1));
        s(i) = str2num(nam(ss(1)+4:ts(1)-1));
    end
end


%%
urs = unique(r);
ucs = unique(c);
uws = unique(w);

sortWs = sort(uws,'ascend');
for cw = 1:length(sortWs)
    wid = sortWs(cw)'
    mif.w(cw).id = wid;  
    allSections = s(w==wid);
    uss = sort(unique(allSections));

    mif.w(cw).sections = uss;
    maxr = max(r(w==wid));
    maxc = max(c(w==wid));
    for cs = 1:length(uss)
        sid = mif.w(cw).sections(cs);
        getSec = (w == wid) & (s == sid);
        mif.w(cw).sec(cs).tileNams = allNams(getSec);
        mif.w(cw).sec(cs).tileFolders = allFolders(getSec);
        mif.w(cw).sec(cs).rowIds = r(getSec);
        mif.w(cw).sec(cs).colIds = c(getSec);
        mosInd = sub2ind([maxr maxc], r(getSec),c(getSec));
        idMos = zeros(maxr,maxc);
        idMos(mosInd) = mosInd;
        mif.w(cw).sec(cs).idMos = idMos;
    end
end
        
%% Identify section files
for w = 1:length(mif.w)
    for s = 1:length(mif.w(w).sec)
       sFold = mif.w(w).sec(s).tileFolders{1};
       mif.w(w).sec(s).secFold = sFold;
       secDir = dir(sFold);
         for d = 1:length(secDir)
             nam = secDir(d).name;
             if ~isempty(regexp(nam,'StageStitched') )
                 mif.w(w).sec(s).ov = nam;
             end
         end
    end
end

mif.Info = imfinfo([allFolders{1} '\' allNams{1}]);


save([TPN 'mif.mat'],'mif')
