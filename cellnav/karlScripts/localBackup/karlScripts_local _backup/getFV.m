function fv = getFV(cidList,fvDir)
if isempty(fvDir)
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
end
fv={};
for i=1:length(cidList)
    curFV=load([fvDir num2str(cidList(i)) '.mat']);
    curFV=curFV.fv;
    curFV.vertices=curFV.vertices(:,[2 3 1]);
    fv{i}=curFV;
end

end