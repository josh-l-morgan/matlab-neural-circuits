SPN = 'E:\IxQ_KarlsRetinaVG3_2019\affFull_fullRes_mipd_single\'

dSPN = dir([SPN '*.png'])
inam = {dSPN.name};

'src'
s = zeros(length(inam),1);
r = s; c = s;
for i = 1:length(inam)
    
    a = sscanf(inam{i},'%d_%d_%d.png')
    s(i) = a(1);
    r(i) = a(2);
    c(i) = a(3);
    
end

maxS = max(s);
maxR = max(r);
maxC = max(c);

fullTileNum  = maxS * maxR * maxC;
tileNum = length(inam);

tileNum/fullTileNum