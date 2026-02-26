function [DotSelect] = axDotSelect(TPN,Dots,in,zum)

DotSelect = Dots;
delIndex = find(in==0);
DotSelect.Pos(delIndex,:) = [];
DotSelect.Vox(delIndex) = [];
DotSelect.Vol(delIndex) = [];
DotSelect.ITMax(delIndex) = [];
DotSelect.ItSum(delIndex) = [];
DotSelect.MeanBrightPuncta(delIndex) = [];
DotSelect.MeanBrightNeurite(delIndex) = [];
DotSelect.Contrast(delIndex) = [];
DotSelect.Num = length(DotSelect.Vox);
DotSelect.DF(delIndex) = [];
DotSelect.DFOf(delIndex) = [];
DotSelect.DFOfTopHalf(delIndex) = [];

clear delIndex
% clearing duplicates
linearPos = sub2ind(DotSelect.ImSize,...
    DotSelect.Pos(:,1),DotSelect.Pos(:,2), DotSelect.Pos(:,3));
for i=1:DotSelect.Num-1
    for j=i+1:DotSelect.Num
        if linearPos(j) == linearPos(i)
            if ~exist('delIndex', 'var')
                delIndex = j;
            else
                next = numel(delIndex)+1;
                delIndex(next) = j;
            end
        else
        end
    end
end

DotSelect.Pos(delIndex,:) = [];
DotSelect.Vox(delIndex) = [];
DotSelect.Vol(delIndex) = [];
DotSelect.ITMax(delIndex) = [];
DotSelect.ItSum(delIndex) = [];
DotSelect.MeanBrightPuncta(delIndex) = [];
DotSelect.MeanBrightNeurite(delIndex) = [];
DotSelect.Contrast(delIndex) = [];
DotSelect.Num = length(DotSelect.Vox);
DotSelect.DF(delIndex) = [];
DotSelect.DFOf(delIndex) = [];
DotSelect.DFOfTopHalf(delIndex) = [];
DotSelect.Strat = std(DotSelect.Pos(:,3)*zum); 

save([TPN 'DotSelect'], 'DotSelect')