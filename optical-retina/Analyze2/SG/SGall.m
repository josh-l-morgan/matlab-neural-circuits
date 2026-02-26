
load(['UseCells.mat'])

for k = 1:size(UseCells,1)
    TPN=char(UseCells(k));
    load([TPN 'find\SG.mat'])
    DT(k)=SG.FinalCrit.DeltaThresh;
    CT(k)=SG.FinalCrit.ConThresh;
    DCT(k)=SG.FinalCrit.DConThresh;
end