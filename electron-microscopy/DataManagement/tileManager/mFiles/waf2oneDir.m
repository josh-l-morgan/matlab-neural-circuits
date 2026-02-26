

SPN = GetMyDir
TPN = [SPN(1:end-1) 'SingleDir\'];
if ~exist(TPN,'dir')
    mkdir(TPN)
end

dSPN = dir(SPN);

for w = 1:length(dSPN)
    
   wafDir = dSPN(w).name;
   if sum(regexp(wafDir,'W'));
    dWPN = dir([SPN wafDir]);
    for s = 1:length(dWPN)
      iNam = dWPN(s).name;
        if sum(regexp(iNam,'.tif'))
           oldName = [SPN wafDir '\' iNam];
               newName = [TPN wafDir '_s' iNam];
               if ~exist(newName,'file')
               for r = 1:3
                   pass = 1;
                   try
                       copyfile(oldName,newName)
                   catch err
                       pass = 0;
                   end
                   if pass 
                       break
                   end
               end
                       
               end
            
        end
        
    end
    
   end
    
    
end
