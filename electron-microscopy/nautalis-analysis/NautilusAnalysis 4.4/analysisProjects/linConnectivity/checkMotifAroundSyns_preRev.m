function[check] = checkMotifAroundSyns(sm,checkS,maxDist)

check.checkS = checkS;
check.maxDist = maxDist;

for i = 1 : length(checkS)
    
    checkPos = sm.pos(checkS(i),:);
    preID = sm.pre(checkS(i));
    postID = sm.post(checkS(i));
    preIDClass = sm.preClass(checkS(i));
    postIDClass = sm.postClass(checkS(i));
    
    check.pos(i,:) = checkPos;
    
    close = find(sm.dist(checkS(i),:)<=maxDist);
    synNum(i) = length(close);
    
    preClass = sm.preClass(close);
    postClass = sm.postClass(close);
    pre = sm.pre(close);
    post = sm.post(close);
    
    %identify all pre and post
    pre2pre = pre(post==preID);
    pre2preClass = preClass(post==preID);
    post2pre = post(pre == preID);
    post2preClass = postClass(pre == preID);
    pre2post = pre(post == postID);
    pre2postClass = preClass(post == postID);
    post2post = post(pre == postID);
    post2postClass = postClass(pre == postID);
    
    %% identify triads
    [recipPre ia ib] = intersect(pre2pre,post2pre);
    [recipPost ia ib] = intersect(pre2post,post2post);
    [io ioa iob] = intersect(pre2pre,post2post);
    ioClass = pre2preClass(ioa);
    [ii iia iib] = intersect(pre2pre,pre2post);
    iiClass = pre2preClass(iia);
    [oo ooa oob] = intersect(post2pre,post2post);
    ooClass = post2preClass(ooa);
    [oi oia oib] = intersect(post2pre,pre2post);
    oiClass = post2preClass(oia);
    
    m.preClass = preIDClass;
    m.postClass = postIDClass;
    
    m.preTcTri = 0;
    m.preLinOutTri =0;
    m.preOnly = 0;
    m.preTcDi = 0;
    m.preLinDi = 0;
    m.preLinInTri = 0;
    
    m.postOnly = 0;
    
    m.reciprocal = 0;
    m.autapse = 0;
    m.recTcTri = 0;
    m.recLinOutTri = 0;
    m.recRgcTri = 0;
    m.recLinInTri = 0;
    m.recLinInOut = 0;
    m.recLinOutIn = 0;
    
    if sum(checkPos>1); %parse motif if in bounds

    %%check pre motif by mapping against post2post
    if isempty(post2post) %local dead end
        m.preOnly = 1; % target has no other output
    elseif preID % not zero
       
        zeroClass = post2postClass(post2post == 0);
        m.preTcDi = m.preTcDi + sum(zeroClass==2);
        m.preLinDi = m.preLinDi + sum(zeroClass==3);
        checkPost = find(post2post>0);
        for o = 1:length(checkPost) % check against all targets of post
            checkPostClass = post2postClass(checkPost(o));
            checkPostID = post2post(checkPost(o));
            if sum(oo == checkPostID); %if preRGC and post LIN converge on target
                if checkPostClass == 2
                    m.preTcTri = m.preTcTri + 1;
                elseif checkPostClass == 3
                    m.preLinOutTri = m.preLinOutTri + 1;
                end
            else % no triad
                if checkPostClass == 2
                    m.preTcDi =  m.preTcDi  + 1;
                elseif checkPostClass == 3
                    m.preLinDi = m.preLinDi  + 1;
                end
            end
        end
    end
    
    %%Find cell that both receives from pre and innervates post
    otherMid = setdiff(oi,0);
    m.preLinInTri = sum(otherMid);
    
    %%determine if post is isolated (otherwise motif will be described from pre perspective)
    if isempty(pre2pre)
        m.postOnly = 1;
    end
    
    %%Check for reciprocity
    m.reciprocal = sum(pre2pre == postID);
    if preID>0
        m.autapse = sum(pre == post);
    end
    if m.reciprocal
        
        %%Recip at top of triad
        shared = find(oo);
        for s = 1:length(shared)
            if ooClass(shared(s)) == 2
                m.recTcTri = m.recTcTri + 1;
            elseif ooClass(shared(s)) == 3;
                m.recLinOutTri = m.recLinOutTri + 1;
            end
        end
        
        %% recip at bottom of triad
        shared = find(ii);
        for s = 1:length(shared)
            if iiClass(shared(s)) == 1
                m.recRgcTri = m.recRgcTri + 1;
            elseif iiClass(shared(s)) == 3;
                m.recLinInTri = m.recLinInTri + 1;
            end
        end
        
        %%Different directions of reciprocal to lin
        shared = find(io);
        m.recLinInOut = length(shared);
        shared = find(oi);
        m.recLinOutIn = length(shared);
        
    end
    end
    check.m(i) = m;
    
end



