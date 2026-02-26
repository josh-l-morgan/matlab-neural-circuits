function[] = writeObjVol(objVol,TPN)


TPNi = [TPN 'objVol\'];


for i = 1:size(objVol,3)
    fileName = sprintf('objVol_%05.0f.png');
    imwrite(objVol(:,:,i),[TPNi fileName])
end


