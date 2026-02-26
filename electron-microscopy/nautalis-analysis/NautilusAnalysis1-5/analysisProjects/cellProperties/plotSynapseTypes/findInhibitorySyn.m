linList = obI.nameProps.cellNum(obI.nameProps.lin)


    preNames = [1028 1032 2033 2034 2035 2003 2007]

    useList.postList
    
    
    for i = 1:length(preNames)
        
        isPre = [synProp.pre] == preNames(i);
        post2pre = [synProp(isPre).post];
        clear isPost
        for t = 1:length(post2pre);
            isPost(t) = sum(linList == post2pre(t))>0
        end
            postLin = post2pre(isPost);
            postNum(i) = sum(isPost);
        
    end
        
        %% 
        
        [synProp([synProp.axoAxonic]>0).pre]
       [synProp([synProp.axoAxonic]>0).post]
