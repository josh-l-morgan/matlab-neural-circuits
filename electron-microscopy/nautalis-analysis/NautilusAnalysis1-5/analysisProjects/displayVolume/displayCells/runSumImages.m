function[] = cells2Dir(obI,dsObj,cellList,Dim,sumDir)


if ~exist(sumDir),mkdir(sumDir); end

%%
fsize = findObSize(dsObj);

for i = 1:length(obI.cell.obIDs)
    cellName = obI.cell.name(i);
    I_Sum = showCellSum(obI,dsObj,cellName,[1,1,1],Dim,fsize)*30;
    image(I_Sum),pause(.01)
    imwrite(I_Sum,[sumDir sprintf('%d.png',cellName)])
end