%% review section List

%load('..\matFiles\us.mat')

%%
noPath = cell(length(us.sec),1);
badTileQual = cell(length(us.sec),1);
meanQual = zeros(length(us.sec),1);
checkByEye = zeros(length(us.sec),1);
usePath = noPath;
parfor s = 1:length(us.sec);
    empties = 0;
    estring = [];
    qualString = [];
    badQuals = 0;
    
    for r = 1:4
        for c = 1:4
            if isempty(us.sec(s).paths{r,c})
                empties = empties+1;
                estring =[estring sprintf('(%d,%d) ',r,c)];
            elseif isempty(usePath{s})
               usePath{s} = us.sec(s).paths{r,c}; 
            end
            
            if us.sec(s).eval.bestQuals(r,c) < 100
               qualString = [qualString sprintf('(%d,%d) ',r,c)];
               badQuals = badQuals+1; 
               usePath{s} = us.sec(s).paths{r,c};
            end
        end
    end
    
    
    
    if empties == 0
        noPath{s} = [];
    elseif empties == 16;
        noPath{s} = 'all';
        checkByEye(s) = 1;

    else
        noPath{s} = estring;
        checkByEye(s) = 1;

    end
    
    if badQuals == 0
        badTileQual{s} = [];

    elseif badQuals == 16;
        badTileQual{s} = 'all';
        checkByEye(s) = 1;
    else
        badTileQual{s} = qualString;
        checkByEye(s) = 1;
    end
    
    meanQual(s) =mean(us.sec(s).eval.bestQuals(:));

    
    if us.sec(s).eval.checkByEye;
        checkByEye(s) = 1;
    end
    if meanQual(s) < 200
        checkByEye(s) = 1;
    end
    
    
    
end


%%

need2check = find(checkByEye);
rev = cell(sum(checkByEye),5);

for s = 1:length(need2check)
    sec = need2check(s);
   
       
       rev{s,1} = [ sec us.sec(sec).wafNum us.sec(sec).wafSec];
       rev{s,2} = usePath{sec};
       rev{s,3} = noPath{sec};
       rev{s,4} = meanQual(sec);
       rev{s,5} = badTileQual{sec};
      
       
end


