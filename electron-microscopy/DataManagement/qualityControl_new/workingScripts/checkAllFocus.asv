function[wQual] = checkAllFocus(WPN)

WPN = GetMyDir
if exist([WPN 'wif.mat'])
    load([WPN 'wif.mat'])
else
wif = getWif(WPN)
end

for s = 1:length(wif.sec)
    nam = [WPN wif.sec(s).name];
    if exist([WPN wif.sec(s).name 'qual.mat'])
        load([WPN wif.sec(s).name 'qual.mat']);
        sec(s) = qual;
    else
        sec(s) = checkSectionQual(nam);
    end
end

%% Build quality image

for s = 1:length(sec)
    mos{s} = sec(s).mos.qualMos;
%     globMos = sec(s).mos.globMos;
%     satMos = sec(s).mos.satMos;
%     mos3{s} = cat(3,qualMos, globMos,satMos);
%     mos3scale{s} = uint8(cat(3, qualMos * 40, globMos * 5, satMos * 20));
     image((mos{s}-1)*100),pause(.1)
end

wQual.mos = mos;
wQual.sec = sec;

save([WPN 'wQual.mat'],'wQual')
