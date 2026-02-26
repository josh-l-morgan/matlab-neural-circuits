TPN = GetMyDir

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
for i = 1:length(matList)
    nam = matList{i};
    rs = regexp(nam,'_r');
    us= regexp(nam,'_');
    cs = regexp(nam,'-c');
    r(i) = str2num(nam(rs(1)+2:cs(1)-1));
    c(i) = str2num(nam(cs(1)+2:us(2)-1));
    load([TPN nam])
    Infos(i) = Info;
    wd(i) = Info.WorkingDistance;
    mos(r(i),c(i)) = Info.WorkingDistance;
end

surf(mos)
mesh(mos)
subplot(1,2,1)
image(fitH(mos))
subplot(1,2,2)
newmos = mos-mean(mos(:));
internalSpread = median(abs(mos(:)-mean(mos(:))));
image(newmos/internalSpread * 50 + 128);
mesh(newmos)
axis([1 7 1 10 internalSpread *-5 internalSpread *5])




