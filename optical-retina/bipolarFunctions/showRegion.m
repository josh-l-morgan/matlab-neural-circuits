%% Get points
% if exist('.\LastRegion.mat')

TPN = GetMyDir;

%%Get resolution
if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
        if isfield(Settings,'ImInfo')
            yxum=Settings.ImInfo{4};
        end
end
if ~exist('yxum'), yxum = 0.0695; end

load([TPN 'Use.mat'])            
            
            
TPNd=dir(TPN); TPNd=TPNd(3:length(TPNd));
Name={};
for i = 1: length(TPNd)
    nam=TPNd(i).name;
    LN = length(nam);
    if LN >3
    if nam(LN-3:LN) == '.rgn'
        Name{length(Name)+1}=nam;
    end
    end
end

TX = load([TPN Name{1}]);
Point=[mean([TX(:,17),TX(:,23)],2) mean([TX(:,18),TX(:,24)],2)];
Point = Point * yxum; %scale for resolution
gains= find(TX(:,4)==65280);  %find all gains (green arrows)
losses = find(TX(:,4) == 255); %find all losses (red arrows)
stable = find(TX(:,4) ==16776960) ; %find all non gain loss puncta
