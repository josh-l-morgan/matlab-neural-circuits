TPN=GetMyDir;

load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])

DrawP=find(SG.passF);

numP=size(DrawP,1);
dp3D=zeros(Dots.ImSize,'uint8');

for i = 1:numP
   
    id=DrawP(i);
    ids=Dots.Vox(id).Ind;
    col=fix(rand*255)+1;
    dp3D(ids)=col;
    
end


dpMax=max(dp3D,[],3);
image(dpMax)