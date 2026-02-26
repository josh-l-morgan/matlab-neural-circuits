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

motDist = 3;

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
    
    
    %% Find Mins
    
    rgcToAll = sm.syn2skel.dist(rgcInTestID,:);
    tcToAll = sm.syn2skel.dist(:,tcInTestID);
    rgcToTC = rgcToAll(:,tcInTestID);
    rgcToMin = mean(min(rgcToTC,[],2)<=motDist);
    tcToMin = mean(min(rgcToTC,[],1)<=motDist);
    
    
    
    reps = 1000;
    skelSize = size(sm.syn2skel.syn2skelLinDist,2);
    numTC = length(tcInTestID);
    rgcToRandMin = zeros(reps,1);
    tcToRandMin = zeros(reps,1);
    rgcToSkel = sm.syn2skel.syn2skelLinDist(rgcInTestID,:);
    for r = 1:reps
        
        randID = ceil(rand(numTC,1)*skelSize);
        rgcToRandTC = rgcToSkel(:,randID);
        rgcToRandMin(r) = mean(min(rgcToRandTC,[],2)<=motDist);
        tcToRandMin(r) = mean(min(rgcToRandTC,[],1)<=motDist);
        
    end
    
    
    
    rgcToRandMin = sort(rgcToRandMin);
    tcToRandMin = sort(tcToRandMin);
    
    rgcToMin
    rgcTo95 = [rgcToRandMin(round(reps*.025)) rgcToRandMin(round(reps * .975))]
    mean(rgcToRandMin)
    P = mean(rgcToRandMin >= rgcToMin)
    
    tcToMin
    tcTo95 = [tcToRandMin(round(reps*.025)) tcToRandMin(round(reps * .975))]
    mean(tcToRandMin)
    P = mean(tcToRandMin >= tcToMin)

    
    
    
    
    
    
    
    
    