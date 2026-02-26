if ~exist('sm','var')
load('D:\LGNs1\Analysis\sm.mat') %% load sm with all relevant synapse data
end

%%
nodeProcType = sm.isTarg + sm.isAx * 2 + sm.isShaft * 3 + sm.isBody * 4;
synProcType = nodeProcType(sm.syn2skel.closestSkel)';

rgcInID = find(sm.syn2skel.preClass == 1);
tcOutID = find(sm.syn2skel.postClass == 2);
rgcInNode = sm.syn2skel.closestSkel(rgcInID);
tcOutNode =  sm.syn2skel.closestSkel(tcOutID);

motDist = 5;
checkType = [1 3];
bins = [0:motDist:40];
nodeLength = sm.skelProps.nodeLength;
clear binRGCtoTC binRGCtoRGC minRGCtoTC minRGCtoRGC binRGCtoNodes  binRGCtoLength
clear   binTCtoNodes binTCtoTC binTCtoRGC minTCtoTC minTCtoRGC binTCtoLength

rgcInTestID = find((sm.syn2skel.preClass == 1) );
    tcInTestID = find((sm.syn2skel.postClass == 2) );
    rgcInTestNode = sm.syn2skel.closestSkel(rgcInTestID);
    tcInTestNode = sm.syn2skel.closestSkel(tcInTestID);
    
    scatter3(sm.skelPos(:,1),sm.skelPos(:,2),sm.skelPos(:,3),'.','k')
    hold on
    scatter3(sm.skelPos(rgcInTestNode,1),sm.skelPos(rgcInTestNode,2),sm.skelPos(rgcInTestNode,3),'o','g')
    scatter3(sm.skelPos(tcInTestNode,1),sm.skelPos(tcInTestNode,2),sm.skelPos(tcInTestNode,3),'o','b')
    hold off
    pause(.1)
    
    
    
    
    for i = 1:length(rgcInTestNode)
        
        %%
        if 0 % show test
            targ = rgcInTestNode(i);
            targS = rgcInTestID(i);
            
            checkPos = sm.skelPos(targ,:);
            b = 20;
            
            scatter3(sm.skelPos(:,1),sm.skelPos(:,2),sm.skelPos(:,3),'.','k')
            hold on
            
            scatter3(sm.skelPos(rgcInTestNode,1),sm.skelPos(rgcInTestNode,2),sm.skelPos(rgcInTestNode,3),'o','g')
            scatter3(sm.skelPos(tcInTestNode,1),sm.skelPos(tcInTestNode,2),sm.skelPos(tcInTestNode,3),'o','b')
            scatter3(sm.syn2skel.synPos(targS,1),sm.syn2skel.synPos(targS,2),sm.syn2skel.synPos(targS,3),'x','r')
            scatter3(sm.skelPos(targ,1),sm.skelPos(targ,2),sm.skelPos(targ,3),'x','b')
            hold off
            ylim([checkPos(2)-b checkPos(2)+b]); xlim([checkPos(1)-b checkPos(1)+b]); zlim([checkPos(3)-b checkPos(3)+b]);
            hold off
            pause
        end
        
        
      %%
%         allDists = sm.skel2skel.linDist(rgcInTestNode(i),:);
%         allDists = sm.skel2skel.dist(rgcInTestNode(i),:);
        allDists = sm.syn2skel.dist(rgcInTestID(i),:);
        
        for b = 1:length(bins)-1
           binRGCtoLength(i,b) =  ...
               sum(nodeLength((allDists>= bins(b)) & (allDists<bins(b+1)))); 
        end
        
        tcDists = allDists(tcOutID);
        rgcDists = allDists(rgcInID);
        binRGCtoTC(i,:) = sum(tcDists<=motDist);
        binRGCtoRGC(i,:) = sum(rgcDists<=motDist);
        
       
        binRGCtoNodes(i,:) = sum(allDists<=motDist);
        minRGCtoTC(i) = min(tcDists);
        minRGCtoRGC(i) = min(rgcDists);
    end
    binRGCtoRGC(i,1) =  binRGCtoRGC(i,1) - 1;
    
    for i = 1:length(tcInTestNode)
        allDists = sm.skel2skel.linDist(tcInTestNode(i),:);
        
        for b = 1:length(bins)-1
           binTCtoLength(i,b) =  ...
               sum(nodeLength((allDists>= bins(b)) & (allDists<bins(b+1)))); 
        end
        
        tcDists = allDists(tcOutNode);
        rgcDists = allDists(rgcInNode);
        
        binTCtoNodes(i,:) = sum(allDists<=motDist);
        binTCtoTC(i,:) = sum(tcDists<=motDist);
        binTCtoRGC(i,:) = sum(rgcDists<=motDist);
        minTCtoTC(i) = min(tcDists);
        minTCtoRGC(i) = min(rgcDists);
    end
    binTCtoTC(i,1) =  binTCtoTC(i,1) - 1;
    



% %% show results
% tCol  = {'g' 'b' 'r'};
% clf,
% subplot(2,1,1)
% 
% hold on
% for t = 1:2
%     
%     
%     densityRGCtoTC = binRGCtoTC./binRGCtoLength;
%         %densityRGCtoTC = binRGCtoTC ;
% 
%     se = std(densityRGCtoTC)./sqrt(size(densityRGCtoTC,1));
%     
%     meanRGCtoTC = mean(densityRGCtoTC,1);
%     plot(bins(2:end),meanRGCtoTC,tCol,'linewidth',3)
%     %     plot(bins(2:end),densityRGCtoTC(round(.05 * size(densityRGCtoTC,1)),:),tCol,'linewidth',2)
%     %     plot(bins(2:end),densityRGCtoTC(round(.95 * size(densityRGCtoTC,1)),:),tCol,'linewidth',2)
%     plot(bins(2:end),meanRGCtoTC+se,tCol,'linewidth',1)
%     plot(bins(2:end),meanRGCtoTC-se,tCol,'linewidth',1)
%     
%     
%     %     meanRGCtoRGC = mean(binRGCtoRGC./binRGCtoLength,1)* motDist;
%     %     plot(bins(2:end),meanRGCtoRGC,tCol)
% end
% 
% hold off
% 
% subplot(2,1,2)
% hold on
% for t = 1:2
%     
%     
%     densityTCtoRGC = binTCtoRGC./binTCtoLength;
%     
%     meanTCtoRGC = mean(densityTCtoRGC,1) ;
%     se = std(densityTCtoRGC)./sqrt(size(densityTCtoRGC,1));
%     
%     plot(bins(2:end),meanTCtoRGC,tCol,'linewidth',3)
%     hold on
%     plot(bins(2:end),meanTCtoRGC+se,tCol,'linewidth',1)
%     plot(bins(2:end),meanTCtoRGC-se,tCol,'linewidth',1)
%     
%     %     meanTCtoTC = mean(binTCtoTC./binTCtoLength,1) * motDist;
%     %     plot(bins(2:end),meanTCtoTC,tCol)
% end
% 
% hold off

    
%% stats

checkB = round(5/motDist);

hitTCtoRGC1 = sum(binTCtoRGC(:,checkB),2);
hitTCtoRGC2 = sum(binTCtoRGC(:,checkB),2);

TChitRGC1 = sum(hitTCtoRGC1>0)
length(hitTCtoRGC1)
TChitRGC2 = sum(hitTCtoRGC2>0)
length(hitTCtoRGC2)
totHit = TChitRGC1 + TChitRGC2



hitRGCtoTC1 = sum(binRGCtoTC(:,checkB),2);
hitRGCtoTC2 = sum(binRGCtoTC(:,checkB),2);

RGChitTC1 = sum(hitRGCtoTC1>0)
length(hitRGCtoTC1)
RGChitTC2 = sum(hitRGCtoTC2>0)
length(hitRGCtoTC2)
totHit = RGChitTC1 + RGChitTC2

sum(binRGCtoTC>0)


