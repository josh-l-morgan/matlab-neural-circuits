function[mos] = grabFocusWB(TPN)

if ~exist('TPN','var')
TPN = GetMyDir
end


%% Find mat files
dTPN = dir(TPN);
matList = {};
for i = 1:length(dTPN)
    nam = dTPN(i).name
    if regexp(nam,'.mat')
        matList{length(matList)+1} = nam;
    end
end

%% load mat files
ct = 0;
for i = 1:length(matList)
    nam = matList{i};
    rs = regexp(nam,'_r');
    us= regexp(nam,'_');
    cs = regexp(nam,'-c');
        if ~isempty(rs) & ~isempty(cs)
            ct = ct + 1;
            r(ct) = str2num(nam(rs(1)+2:cs(1)-1));
            c(ct) = str2num(nam(cs(1)+2:us(2)-1));
            load([TPN nam])
            Infos(ct) = Info;
            wd(ct) = Info.WorkingDistance;
            mos(r(ct),c(ct)) = Info.WorkingDistance;
        end

end
