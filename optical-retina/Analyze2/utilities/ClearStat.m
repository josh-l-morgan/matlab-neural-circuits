%% Collect all relevant information about cells


%Get directory name
KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));


%Run cells to clear a stat
for k = 1:size(Kdir,1)
    
    TPN = [KPN  Kdir(k).name '\'];
    DPN = [TPN 'I\']   
    DPNd=[TPN 'data\'];
    
    
   
   
   if exist([TPN 'SKStatus.mat'])
       load([TPN 'SKStatus.mat'])
       SKStatus.FM=0;
       save([TPN 'SKStatus.mat'],'SKStatus')

   end  
   clear SKStatus;
   

end




