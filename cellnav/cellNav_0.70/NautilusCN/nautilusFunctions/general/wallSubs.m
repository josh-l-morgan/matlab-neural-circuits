function[newSub] = wallSubs(subs,wallVec);

%%Replace values in subs esceding the wall limits with the wall limits



if size(subs,2) == 3
    
    newSub = subs;
    if size(wallVec,1) == 2;
        newSub(subs(:,1)< wallVec(1,1),1) = wallVec(1,1);
        newSub(subs(:,2)< wallVec(1,2),2) = wallVec(1,2);
        newSub(subs(:,3)< wallVec(1,3),3) = wallVec(1,3);
    end
    
    newSub(subs(:,1)> wallVec(end,1),1) = wallVec(end,1);
    newSub(subs(:,2)> wallVec(end,2),2) = wallVec(end,2);
    newSub(subs(:,3)> wallVec(end,3),3) = wallVec(end,3);
    
else
    newSub = subs;
    %disp('wrong size');
end