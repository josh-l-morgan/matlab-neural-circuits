function[skelOri pc] = princompSkeletons(cellStruc,seedSubRaw,lookRange);

%voxSize = [.512 .512 .480]; % in microns
%
% anchorScale = [.0184 0.016 0.030];
% voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
% vastScale = [1 4/4.6 1];
% pixum = .0184 * 8;
%
% voxSize = voxelScale *2;
arbor = cellStruc.arbor;
voxSize = cellStruc.voxSize;
seedSub = seedSubRaw;

if ~exist('lookRange','var')
    lookRange = [5 10]; % min and max value to look at
end
subplotY = 4;

%%
%seedSub = arbor.seed .* voxSize;
allEdgeMids = [];
allEdgeLengths = [];
for b = 1:length(arbor.branches)
    allEdgeMids = cat(1,allEdgeMids,arbor.branches(b).edgeMids);
    allEdgeLengths = cat(1,allEdgeLengths,arbor.branches(b).edgeLengths');
end

for d = 1:3
    allEdgeMids(:,d) = allEdgeMids(:,d) * voxSize(d);
end
edgePos = cat(2,(allEdgeMids(:,1)-seedSub(:,1)), ...
    (allEdgeMids(:,2)-seedSub(:,2)),(allEdgeMids(:,3)-seedSub(:,3)) );

edgeDist = sqrt(edgePos(:,1).^2 + edgePos(:,2).^2  + edgePos(:,3).^2);
inRange = (edgeDist>lookRange(1)) & (edgeDist<=lookRange(2));


usePos = edgePos(inRange,:);
[Coef,Score,latent,tsquare] = princomp(usePos);
skelOri.prinSTD = std(Score);
skelOri.cartSTD = std(edgePos);

scatter(Score(:,1),Score(:,2),'.')
hold on
%scatter(Transformed(:,1),Transformed(:,2),'.','r')
hold off
xlim([-50 50])
ylim([-50 50])

pc.Score = Score;
pc.Coef = Coef;
pc.latent = latent;
pc.tsquare = tsquare;
%scatter(edgePos(:,1),edgePos(:,3),'.')

%edgePos = Score;

%% cardinals
card.help = ['Analyizes points within look range in terms of cardinal directions. ' ...
    'results come in the form of rations of variance and ratios of points assigned a '...
    'cardinal direction'];
card.var = [var(usePos(:,1)) var(usePos(:,2)) var(usePos(:,3))];
card.varRat = [card.var(2)/card.var(3) card.var(1)/card.var(3)];
is3rd2nd = sum((abs(usePos(:,2))> abs(usePos(:,3))))/sum((abs(usePos(:,3))> abs(usePos(:,2))));
is3rd1st = sum((abs(usePos(:,1))> abs(usePos(:,3))))/sum((abs(usePos(:,3))> abs(usePos(:,1))));
card.upSideFront = [is3rd2nd is3rd1st];



%%


useScore= Score;
useEdgeLengths = allEdgeLengths(inRange);

dirCount = zeros(3,2);
for e = 1:size(useScore,1)
    edg = useScore(e,:);
    dimTarg = find(abs(edg) == max(abs(edg)));
    
    if edg(dimTarg)>0
        dirTarg = 1;
    else
        dirTarg = 2;
    end
    dirCount(dimTarg,dirTarg) = dirCount(dimTarg,dirTarg) + ...
        useEdgeLengths(e);
end

%     subplot(subplotY,3,[7:9])
%     bar(dirCount(:,1)*-1)
%     hold on
%     bar(dirCount(:,2))
%     hold off
%

skelOri.latent = latent;
skelOri.latRat = [(latent(1)-latent(2))/latent(1) (latent(1)-latent(3))/latent(1)];
skelOri.dirCount = dirCount;
ori = sum(dirCount,2);
skelOri.oriCount = ori;
skelOri.oriRat = [(ori(1)-ori(2))/(ori(1))...
    (ori(1)-ori(3))/(ori(1))];
skelOri.lengthUsed = sum(useEdgeLengths);
skelOri.totalLength = sum(allEdgeLengths);
skelOri.lookRange = lookRange;
skelOri.card = card;



%% Classify by model








