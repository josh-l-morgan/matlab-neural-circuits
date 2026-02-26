%{
branchDat = {};
synDat = {};
%}


branchPos = zeros(size(branchDat,1),3);
branchType = zeros(size(branchDat,1),1);
for i = 1:size(branchDat,1)
    pos = branchDat{i,1};
    tag = branchDat{i,2};
    wasNum = 0;
    d = 0;
    tempn = [];
    for s = 1:length(pos) %parse possition
        
        n = str2num(pos(s));
        if ~isempty(n)
            if ~wasNum
                d = d+1; %increment next dimension
                tempn = []; %make new temp
            end
            tempn = [tempn pos(s)];
            
            wasNum = 1;
        else
            if wasNum
                branchPos(i,d) = str2num(tempn);
            end
            wasNum = 0;
        end
    end
    
    if ischar(tag)
        
        if sum(regexp(tag,'thick'))
            branchType(i) = 1;
        elseif sum(regexp(tag,'thin'))
            branchType(i) = 2;
        else
            branchType(i) = 0;
        end
    else
        branchType(i) = 0;
    end
    
    
end



%%

synPos = zeros(size(synDat,1),3);
synType = zeros(size(synDat,1),1);
for i = 1:size(synDat,1)
    pos = synDat{i,1};
    tag = synDat{i,2};
    wasNum = 0;
    d = 0;
    tempn = [];
    for s = 1:length(pos) %parse possition
        
        n = str2num(pos(s));
        if ~isempty(n)
            if ~wasNum
                d = d+1; %increment next dimension
                tempn = []; %make new temp
            end
            tempn = [tempn pos(s)];
            
            wasNum = 1;
        else
            if wasNum
                synPos(i,d) = str2num(tempn);
            end
            wasNum = 0;
        end
    end
    
    if ischar(tag)
        
        
        if sum(regexp(tag,'pre'))
            synType(i) = 1;
        elseif sum(regexp(tag,'pos'))
            synType(i) = 2;
        else
            synType(i) = 0;
        end
    else
        synType(i) = 0;
        
    end
    
    
end





%%
branchList1 = sum(branchPos,2);
branchList2 = sum(branchPos,2) * 0;

showBranch = branchPos(branchList1>0,:);


preList = sum(synPos,2) & (synType == 1);
postList = sum(synPos,2) & (synType == 2);

showPre = synPos(preList>0,:);
showPost = synPos(postList>0,:);

scatter3(showBranch(:,1),showBranch(:,2),showBranch(:,3),'k','.')
hold on
scatter3(showPre(:,1),showPre(:,2),showPre(:,3),'g','.')
scatter3(showPost(:,1),showPost(:,2),showPost(:,3),'r','.')
hold off







