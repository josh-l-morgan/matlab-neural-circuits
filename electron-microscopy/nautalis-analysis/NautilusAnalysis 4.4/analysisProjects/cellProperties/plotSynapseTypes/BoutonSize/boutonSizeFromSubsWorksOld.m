%%Make predictions of how many synapses each pair of cells should form
%%based of the distribution of their skeletons. -
%-Skeletons are filtered



%% Plot axon to tcr
clear all
%MPN = GetMyDir;
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

seedList = [108 201 109 903 907]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
preList = useList.preList;
colormap gray(256)


%%

anchorScale = [.0184 0.016 0.030]*4;
voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3)* 4 ];


%% Get cells

axList = preList;
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);

for r = 1:50;
   ball = ones(r*2,r*2,r*2);
   ballInd = find(ball);
   [y x z] = ind2sub(size(ball),ballInd);
   
   dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+(z-mean(z)).^2);
   ball(dists>r) = 0; 
   ball = ball/sum(ball(:));
   image(squeeze(sum(ball,2))*100000); 
   balls{r} = ball;
end

%%  Get Skeletons
 colormap jet(256)
checkRad = 50;
for a = 1:length(axList)
        disp(sprintf('running %d of %d',a,length(axList)))
            subs = scaleSubs(getCellSubs(obI,dsObj,axList(a)),voxelScale);
    
            isAx = find((edges(:,2) == axList(a)) & sum(synAnchors,2));
            
            scatter(subs(:,1),subs(:,2),'.')
            hold on
            scatter(synAnchors(isAx,1),synAnchors(isAx,2),'r','.')
            hold off
            sc = 0;
            clear butProfile synPos
            for s = 1:length(isAx);
                targSyn = isAx(s);
                preDists = sqrt((subs(:,1)-synAnchors(targSyn,1)).^2 + (subs(:,2)-synAnchors(targSyn,2)).^2 + ...
                    (subs(:,3)-synAnchors(targSyn,3)).^2);
                minDist = min(preDists);
                targPre = find(preDists == minDist,1);
                
                    nearVox = find(preDists<=checkRad);
                    if ~isempty(nearVox)
                    sc = sc + 1;
                    
                    targPre = find(nearVox == targPre,1);
                    nearDist = preDists(nearVox);
                    nearSubs = subs(nearVox,:);
                    minVox = min(nearSubs,[],1);
                    nearSubs = round(nearSubs - repmat(minVox,[size(nearSubs,1) 1]))+1;
                    maxVox = max(nearSubs,[],1);
                    nearVol = zeros(maxVox);
                    nearInd = sub2ind(maxVox,nearSubs(:,1),nearSubs(:,2),nearSubs(:,3));
                    nearVol(nearInd) = 1;
                    image(squeeze(sum(nearVol,1))*10)
                    
                    synPos(sc,:) = synAnchors(targSyn,:);
                    for b = 1:20;
                        filtVol = nearVol;
                        ball = balls{b};
                        filtVol(nearInd(nearDist>(b*2))) = 0;
                        %Ic = convn(filtVol,ball,'same');
                        Ic = fastCon(filtVol,ball);

                        butProfile(sc,b) = max(Ic(:));
                        
%                         maxIc = max(Ic,[],3);
%                         
%                         subplot(3,1,1)
%                         image(sum(nearVol,3)*10);
%                         subplot(3,1,2)
%                         image(sum(ball,3)*256/max(ball(:))/b/2)
%                         ylim([1 size(nearVol,1)])
%                         xlim([1 size(nearVol,2)])
%                         subplot(3,1,3)
%                         image(maxIc*256)
%                         pause(.01)
                    end
                    
                    
                end
            end
            
            axProfile{a} = butProfile;
            mapSyn{a} = synPos;
            
          
end

%% Analyize profiles
fitThresh = 0.5; %define fit to ball

ballVol = zeros(length(balls),1);
for b = 1:length(balls)
    ball = balls{b};
    ballVol(b) = sum(ball(:)>0);
end

histButVolBin = [0:100:3000];
histButVol = zeros(length(axProfile),length(histButVolBin));
for a = 1:length(axProfile)
    a
prof = axProfile{a};
plot(prof'),pause(.01)
    goodProf = 0;
    clear butVol
    for p = 1:size(prof,1);
        checkProf = prof(p,:);
        ballFit = max(find(checkProf> fitThresh));
        
        if ~isempty(ballFit)
            goodProf = goodProf + 1;
            butVol(goodProf) = ballVol(ballFit)*checkProf(ballFit);
            
        end
    end
    butVols{a} = butVol;
    histButVol(a,:) = hist(butVol,histButVolBin);    

end

subplot(2,1,1)
image(histButVol*4)
subplot(2,1,2)











