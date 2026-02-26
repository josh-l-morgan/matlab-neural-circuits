function [DotSelect] = axDotSelect2(TPN,Dots,lowBound,highBound,zum)

DotSelect = Dots;

delIndex = find(Dots.MeanBright./Dots.DFOf > lowBound |...
    Dots.MeanBright./Dots.DFOf < highBound);

DotSelect.Pos(delIndex,:) = [];
DotSelect.Vox(delIndex) = [];
DotSelect.Vol(delIndex) = [];
DotSelect.ITMax(delIndex) = [];
DotSelect.ItSum(delIndex) = [];
DotSelect.MeanBrightPuncta(delIndex) = [];
DotSelect.MeanBrightNeurite(delIndex) = [];
DotSelect.Num = length(DotSelect.Vox);
DotSelect.DF(delIndex) = [];
DotSelect.DFOf(delIndex) = [];
DotSelect.DFOfTopHalf(delIndex) = [];
DotSelect.Strat = std(DotSelect.Pos(:,3)*zum); 

save([TPN 'DotSelect'], 'DotSelect')