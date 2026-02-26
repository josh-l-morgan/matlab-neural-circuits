

load('UseCells')

for k = 1:size(UseCells,1)
    k
   TPN = char(UseCells(k));
   load([TPN 'data\DepthDevNoCB.mat'])
   save([TPN 'data\DepthDevNoCB.mat'],'DepthDev')
    
    
end