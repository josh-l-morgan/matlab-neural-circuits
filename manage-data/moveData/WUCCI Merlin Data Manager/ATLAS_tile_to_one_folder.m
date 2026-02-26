


SPN = GetMyDir
TPN = [SPN(1:end-1) '_tile_r1_c2\'];
if ~exist(TPN), mkdir(TPN),end


fSPN = findFolders(SPN)


%%
clear iFold iNam ids iT
tifNum = 0;
for f = 1:length(fSPN)
    fold = fSPN{f};
    dFold = dir(fold);
    for d = 1:length(dFold)
        nam = dFold(d).name;
        if sum(regexp(nam,'.tif'))
            tifNum = tifNum+1;
            r = regexp(nam,'_r');
            c = regexp(nam,'-c');
            s = regexp(nam,'_S_');
            u = regexp(nam,'_');
            
            R = str2num(nam(r+2:c-1));
            C = str2num(nam(c+2:s-1));
            S = str2num(nam(s+3:u(end)-1));
            T = str2num(nam(u(end)+1:end-4));
            id = S * 1000000000 + R*1000000 + C * 1000 ;
            
            iFold{tifNum} = fold;
            iNam{tifNum} = nam;
            ids(tifNum,:) = [ R C S];
            iT(tifNum) = T;
        end
    end
end
   

%% get unique sections
uS = unique(ids(:,3));
useTimes = []; %record times to use
clear secFiles
for i = 1:length(uS)
   isSec = find(ids(:,3) == uS(i));
   folds = iFold(isSec);
   found{1} = folds{1};
   foldSeg = zeros(1,length(folds));
   numSeg = 0;
   clear segFold
   for f = 1:length(folds)-1;
       if ~foldSeg(f)
            numSeg = numSeg+1;
            foldSeg(f) = numSeg;
            segFold{numSeg} = folds{f}
       end
       nam = folds{f};
       for r = f+1:length(folds);
           if ~foldSeg(r) %if unassigned
               if strcmp(nam,folds{r});
                   foldSeg(r) = foldSeg(f);
               end
           end
       end
       
       
   end
   foldSeg
    
      
         times = iT(isSec) ;
         bestTimes = find(times == max(times));
         useSec = isSec(bestTimes);
    
    secDat(i).secID = uS(i);
    secDat(i).folder = iFold{useSec(1)};
    secDat(i).secFold = sprintf('sec_%06.0f',uS(i));
    secDat(i).iNams = iNam(useSec);
    secDat(i).ids = ids(useSec,:);
    secDat(i).mainDir = TPN;
end


%% copy Data
save([TPN 'secDat.mat'],'secDat')
for s = 1:length(secDat)
    for i = 1:length(secDat(i).iNams)
        
        try    if (secDat(s).ids(i,1) == 1) & (secDat(s).ids(i,2) == 2)
               
                
        disp(sprintf('copying section number %d of %d',s,length(secDat)))
            
            oldName = [secDat(s).folder '\' secDat(s).iNams{i}]
            secFold = [secDat(s).mainDir];
            if ~exist(secFold,'dir'),mkdir(secFold);end
            newName = [secFold secDat(s).iNams{i}];
            if ~exist(newName,'file')
                copyfile(oldName,newName);
            end
            end
        end
        
    end
 end


    
    

