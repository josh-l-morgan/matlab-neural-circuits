function[mif] = getMif(TPN)

if ~exist('TPN','var')
    TPN = GetMyDir;
end

APN = findFolders(TPN);
    allNams = {}; allTimes = []; allFolders = {};
for f = 1: length(APN)
    [aNams aTimes]= getTiles(APN{f});
    fID = ones(length(aNams),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
    allNams = cat(1,allNams, aNams);
    allTimes = cat(1,allTimes, aTimes);
end


wif.dir = TPN;
% 
% dTPN = dir(TPN);
% ct = 0;
% for i = 1:length(dTPN);
%     nam = dTPN(i).name;
%     findSec = regexp(nam,'_Sec');
%     afterSec = regexp(nam,'_Montage');
%     if ~isempty(findSec) & ~isempty(afterSec)
%         ct = ct + 1;
%         wif.sec(ct).num = str2num(nam(findSec+4:afterSec-1));
%         wif.sec(ct).name = [nam '\'];
%     end
% end


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
            r = str2num(nam(rs(1)+2:cs(1)-1));
            c = str2num(nam(cs(1)+2:us(2)-1));
            w = str2num(nam(ws(1)+2:ss(1)-1));
            s = str2num(nam(ss(1)+4:ts(1)-1));
            
            ct = ct+1;
            wif.sec(s).tile(ct).row = r;
            wif.sec(s).tile(ct).col = c;
            wif.sec(s).tile(ct).name = nam;
            load([TPN  secNam nam(1:end-4) '.mat'])
            wif.sec(s).tile(ct).Info = Info;            
        elseif ~isempty(regexp(nam,'StageStitched'))
            wif.sec(s).ov =  nam;
        end
    end
end

save([TPN 'wif.mat'],'wif')
