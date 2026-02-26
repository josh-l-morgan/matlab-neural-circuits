function[check] = getMotifsAroundSyn(sm,checkS,maxDist,targCell)

%%check motifs relative to a list of synapses on sm

check.checkS = checkS;
check.targCell = targCell;
check.maxDist = maxDist;

for i = 1:length(checkS) % Define position by input synapse
    checkPos = sm.pos(checkS(i),:);
    preID = sm.pre(checkS(i));
    postID = sm.post(checkS(i));

    check.pos(i,:) = checkPos;
    
    close = find(sm.dist(i,:)<=maxDist);
    synNum(i) = length(close);
    if isempty(close)
        mot{i} = 0;
    else
        
        preClass = sm.preClass(close);
        postClass = sm.postClass(close);
        pre = sm.pre(close);
        post = sm.post(close);
  
        
%         rgcIn = (preClass == 1) & (post == targCell);
%         tcOut = (postClass == 2) & (pre == targCell);
%         linIn = (preClass == 3) & (post == targCell);
%         linOut = (postClass == 3) & (pre == targCell);
%         unkIn = (preClass == 4) & (post == targCell);
        
        m.preLab = length(checkS);
        m.preTc = 0;
        m.preTcDoub = 0;
        m.preLin = 0;
        m.preLinTc = 0;
        m.preLinLinIn = 0;
        m.preLinLinOut = 0;
        m.preLinLinDoub = 0;
        m.preOneOut = 0;
        m.preLinTcDi = 0;
        m.preLinLinDi = 0; %chain of two
        m.recip = 0;
        m.autapse = 0;
        
        
        %% check from Pre perspective
        isID = pre==preID;
        rPostID = post(isID);
        rPostClass = postClass(isID);
        
        if isempty(setdiff(rPostID,postID))
            m.preOneOut = m.preOneOut + 1;
        end
                
        %%check if pre is reciprical or autapse
        m.recip = m.recip + sum((pre==targCell) & (post==preID));
        if preID == targCell
            m.autapse = m.recip;
        end
        
        m.preTc = m.preTc + sum(rPostClass == 2);
        m.preLin = m.preLin + sum(rPostClass == 3);
        
        %Add zeros
        m.preTcDoub = m.preTcDoub + sum((rPostClass==2) & (rPostID==0));
        m.preLinLinDoub = m.preLinLinDoub + sum((rPostClass==3) & (rPostID==0));
        
        %check triad vs double of labeled cells
        checkPost = setdiff(unique(rPostID),[0 postID]);
        
        %%find lin triads
        for p = 1:length(checkPost)
            
            if sum((post == checkPost(p)) & (postClass == 2)) % if post is TC
                isPost = (pre==postID) & (post == checkPost(p)) ; %where TC also is innervated by target
                if sum(isPost)
                    m.preLinTc = m.preLinTc + 1;
                else
                    m.preTcDoub = m.preTcDoub + 1;
                end
            end
            
            if sum((post == checkPost(p)) & (postClass == 3)) %if post is LIN
                isPre = (pre == checkPost(p)) & (post == postID)& (preClass == 3); %where lin innervates target
                if sum(isPre)
                    m.preLinLinIn = m.preLinLinIn + 1;
                else
                    m.preLinLinDoub = m.preLinLinDoub +1;
                end
            end
            
            if sum((post == checkPost(p)) & (postClass == 3))%if post is LIN
                isPost = (pre==postID) & (post == checkPost(p)) & (postClass == 3); %where lin is innervated by target
                if sum(isPost)
                    m.preLinLinOut = m.preLinTc + 1;
                else
                    m.preLinLinDoub = m.preLinLinDoub +1;
                end
            end
            
        end
        
        %%Find RGC LIN TC diads
        postTarg = post((pre == postID) & (postClass == 2));
        m.preLinTcDi = m.preLinTcDi + sum(postTarg == 0);
        m.preLinTcDi = m.preLinTcDi + length(setdiff(postTarg,[0; checkPost]));
        
        %%Find RGC LIN LIN diads
        postTarg = post((pre == postID) & (postClass == 3));
        m.preLinLinDi = m.preLinLinDi + sum(postTarg == 0);
        m.preLinLinDi = m.preLinLinDi + length(setdiff(postTarg,[0; checkPost]));
        
        %% check from post perspective
        
        
        
        
        
        
        
    end
    check.m(i) = m;
    
end

