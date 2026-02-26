



SPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2_reDice7\'
TPN = 'X:\Active\morganLab\karlsRetina\scm_ixQ_hires_pyr2_reDice7\'

L = length(SPN);

while 1
           
   files = dir([SPN '**\*'])
   
   for f = 1:length(files)
        
       
       oldFile = [files(f).folder '\' files(f).name];
       newDir = [TPN files(f).folder(L+1:end) '\' ];
       newFile = [newDir files(f).name];
       
       
        if ~exist([files(f).folder files(f).name])
            if ~exist(newDir,'dir')
                mkdir(newDir);
            end
            try
                copyfile(oldFile,newFile)
            end
        end
              
   end
   pause(10)
           
end
    
    





