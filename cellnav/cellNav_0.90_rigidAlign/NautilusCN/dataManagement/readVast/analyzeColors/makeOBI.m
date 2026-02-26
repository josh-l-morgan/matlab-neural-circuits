function[obI] = makeOBI(TPN,datSheet)

load('MPN.mat')

obI.colStruc = readVastColors2019;
obI = fixColStruc(obI);

if exist([MPN 'dat.mat'],'file')
    load([MPN 'dat.mat'])
    obI.dat = parseGoogleDat;
end

obI.nameProps = getNameProps2019(obI.colStruc.names);
obI = getObiCellProps(obI);



save([TPN 'obI.mat'],'obI');














