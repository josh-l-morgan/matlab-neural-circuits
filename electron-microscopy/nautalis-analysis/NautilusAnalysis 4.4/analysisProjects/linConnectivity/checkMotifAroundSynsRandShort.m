function[check] = checkMotifAroundSynsRandShort(sm,checkS,maxDist)

check.checkS = checkS;
check.maxDist = maxDist;

% scatter3(sm.skelPos(:,1),sm.skelPos(:,2),sm.skelPos(:,3),'.','c')
%     hold on
for i = 1 : length(checkS)
    
    
    checkPos = sm.pos(checkS(i),:);
    preID = sm.pre(checkS(i));
    postID = sm.post(checkS(i));
    preIDClass = sm.preClass(checkS(i));
    postIDClass = sm.postClass(checkS(i));
    
    check.pos(i,:) = checkPos;
    
    
    %%Make Random Distances
    
    distToSkel = sm.skelDist(:,checkS(i));
    distR = distToSkel(ceil(rand(length(sm.pre),1)*length(distToSkel))); %pick a random skeleton position for every synapse
   
    close = find(distR<=maxDist);
    %close = find(sm.dist(checkS(i),:)<=maxDist); %Not random?
    synNum(i) = length(close);
    
    synPos = sm.pos(close,:);
%     scatter3(checkPos(1),checkPos(2),checkPos(3),10,'k');
%     scatter3(synPos(:,1),synPos(:,2),synPos(:,3),'r');
%     xlim([checkPos(1)-10 checkPos(1)+10])
%     ylim([checkPos(2)-10 checkPos(2)+10])
%     zlim([checkPos(3)-10 checkPos(3)+10])
% 
%     pause(.01)
    
    preClass = sm.preClass(close);
    postClass = sm.postClass(close);
    pre = sm.pre(close);
    post = sm.post(close);
    
    %check without identity
    m.closePostTC = sum((postClass==2));
    m.closePreRGC = sum((preClass==1));
    m.closePostLIN = sum((postClass==3));
    m.closePreLIN = sum((preClass==3));
    
   
    m.diadTCs = sum((postClass==2) & (pre == postID)); %count target to TC outputs
    m.diadRGCs = sum((preClass==1) & (post == preID)); %count RGC to presynaptic cell inputs 
    
    check.m(i) = m;
    
end



