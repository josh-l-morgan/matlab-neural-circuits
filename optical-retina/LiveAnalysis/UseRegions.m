clear all


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


%% relative distance to gain or loss

%%find distance to nearest gain

for i = 1:size(Point,1)
    nG = dist(Point(gains,:),Point(i,:));
    nearestGain(i) = min(nG(nG>0));    
end
nBins = 0:2:max(nearestGain)
hist(nearestGain(stable),nBins,'b')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')
hold on
hist(nearestGain(gains),nBins,'FaceColor','r')
hold off



if length(losses)>1
    for i = length(stable)
    nL = dist(Point(losses,:),Point(i,:));
    nearestLoss(i) = min(nL(nL>0));
    end
end

%% Probability recovery profile

numDot = size(Point,1);
distMat = zeros(numDot,numDot);
for i = 1:numDot
    distMat(i,:) = dist(Point,Point(i,:));
end


lookBack = 100;
%plot(1),hold on
for i = 1: size(Point,1)
     distMids = dist([Use.Mids(:,1) Use.Mids(:,2)],Point(i,:));
    for d = 1:100
        nearLength = sum(Use.Length(distMids<=d));
        nearStable = sum((distMat(i,stable)<=d) & sum(distMat(i,stable)<=d));
        nearGain = sum(distMat(i,gains)<=d) & sum(distMat(i,stable)<=d);
        nearLoss = sum(distMat(i,losses)<=d) & sum(distMat(i,stable)<=d);
        pRecov(i,d) =  nearGain/nearLength;
    end
%     if ~isempty(find(gains==i))
%         plot(pRecov(i,:),'g')
%     elseif ~isempty(find(stable==i))
%         plot(pRecov(i,:),'b')
%     else
%         plot(pRecov(i,:),'r')
%     end
%     pause
    
end
hold off


plot(mean(pRecov(gains,:),1),'g')
hold on
plot(mean(pRecov(stable,:),1),'b')
plot(mean(pRecov(losses,:),1),'r')
hold off

plot(mean(pRecov(gains,:),1)-mean(pRecov(stable,:),1))


%% Poisson
%pRecov = pRecov(:,6:size(pRecov,2));
rate = mean(pRecov,1);
num = mean(pRecov(gains,:),1);
(rate.^num .* exp(-rate))./factorial(num);


