function[col] = showSkel(skel);


%%
subs = skel.nodeSubs;
for i = 1:3
    subs(:,i) = subs(:,i) - min(subs(:,i)) + 1;
end
maxSubs = max(subs,[],1)+1;

predIm  = zeros(maxSubs(1),maxSubs(2));

usedTip = zeros(length(skel.nodes),1);
for t = 1:length(skel.tips)
    tip = skel.tips(t);
    for r = 1:length(skel.nodes)
        %if usedTip(tip), break,end %stop when you get to a previously used vox
        tSub = subs(tip,:);
        predIm(tSub(1),tSub(2)) = 1;% predIm(tSub(1),tSub(2)) + 1;
        tip = skel.pred(tip);
        if tip<1,break,end
        usedTip(tip) = 1;
    end
end

image(predIm*100)


%% draw edges
edgeIm = predIm * 0;
edges = skel.edges;
drawNum = 10;
steps = [0:drawNum]/drawNum;
stepE = zeros(length(steps),3);

for i = 1:size(edges,1);
   startE = subs(edges(i,1),:);
   stopE = subs(edges(i,2),:);
   difE = [(stopE(1)-startE(1))  (stopE(2)-startE(2)) ...
       (stopE(3)-startE(3))];
   %lengthE = sqrt(sum(difE.^2));
   
   stepE(:,1) = startE(1) + difE(1) * steps ;
   stepE(:,2) = startE(2) + difE(2) * steps ;
   stepE(:,3) = startE(3) + difE(3) * steps ;
   stepE = round(stepE);
   
   indE = sub2ind(size(edgeIm),stepE(:,1),stepE(:,2));
   for s = 1:length(indE)
       
       edgeIm(indE(s)) = edgeIm(indE(s)) + 1;
   end
end

image(edgeIm*100)

%% Draw bridges

%% draw edges
bridgeIm = predIm * 0;
edges = skel.bridge;
drawNum = 100;
steps = [0:drawNum]/drawNum;
stepE = zeros(length(steps),3);
for i = 1:size(edges,1);
   startE = subs(edges(i,1),:);
   stopE = subs(edges(i,2),:);
   difE = [(stopE(1)-startE(1))  (stopE(2)-startE(2)) ...
       (stopE(3)-startE(3))];
   %lengthE = sqrt(sum(difE.^2));
   
   stepE(:,1) = startE(1) + difE(1) * steps ;
   stepE(:,2) = startE(2) + difE(2) * steps ;
   stepE(:,3) = startE(3) + difE(3) * steps ;
   stepE = round(stepE);
   
   indE = sub2ind(size(edgeIm),stepE(:,1),stepE(:,2));
   for s = 1:length(indE)
    bridgeIm(indE(s)) = bridgeIm(indE(s)) + 1;
   end
    
end
    
image(bridgeIm*100)

%%

col = bridgeIm * 100;
col(:,:,2) = predIm*300;
col(:,:,3) = edgeIm * 150/median(edgeIm(edgeIm>0));

image(uint8(col)),pause(.010)
 


