



SPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2_reDice3\'

L = length(SPN);

           
   files = dir([SPN '**\*'])
   
   for f = 1:length(files)
       if ~mod(f,1000), disp(sprintf('%d of %d',f,length(files)));end
       if ~files(f).isdir
       delFile = [files(f).folder '\' files(f).name];  
            try
                delete(delFile)
            end
            
       end
   end
   
    files = dir([SPN '**\*'])

   for f = 1:length(files)
       disp(sprintf('%d of %d',f,length(files)))
       if files(f).isdir
       delFold = [files(f).folder '\' files(f).name];  
            try
                rmdir(delFold)
            end
            
       end
   end
   
   
   rmdir(SPN)
    





