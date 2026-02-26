%% Collect all relevant information about cells


%Get directory name
KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));

%Create Cell Arrays
S=10;
SS=cell(size(Kdir,1),S);
SSheader={'Name' 'DF' 'Ra' 'Ma' 'Rd' 'Crit' 'FM' 'AllSeg' 'AllSegCut'}

SOut=cell(size(Kdir,1),4);
SOheader={'Name' 'Age' 'Type' 'AverageDensity'};

%Run cells to gather info
for k = 1:size(Kdir,1)
    
    TPN = [KPN  Kdir(k).name '\'];
    DPN = [TPN 'I\']   
    DPNd=[TPN 'data\'];
    
    
    %Cell Name
    SS(k,1)=cellstr(Kdir(k).name);
   
   if exist([TPN 'SKStatus.mat'])
       load([TPN 'SKStatus.mat'])
       if isfield(SKStatus,'DF')
           SS(k,2)=num2cell(SKStatus.DF,2);
       else
           SS(k,2)=num2cell('NA',2); 
       end       
       if isfield(SKStatus,'Ra')
           SS(k,3)=num2cell(SKStatus.Ra,2);
       else
           SS(k,3)=num2cell('NA',2); 
       end 
       if isfield(SKStatus,'Ma')
           SS(k,4)=num2cell(SKStatus.Ma,2);
       else
           SS(k,4)=num2cell('NA',2); 
       end
       if isfield(SKStatus,'Rd')
           SS(k,5)=num2cell(SKStatus.Rd,2);
       else
           SS(k,5)=num2cell('NA',2); 
       end
       if isfield(SKStatus,'Crit')
           SS(k,6)=num2cell(SKStatus.Crit,2);
       else
           SS(k,6)=num2cell('NA',2); 
       end
       if isfield(SKStatus,'FM')
           SS(k,7)=num2cell(SKStatus.FM,2);
       else
           SS(k,7)=num2cell('NA',2); 
       end
   end  
   clear SKStatus
   
   SS(k,8)= num2cell(exist([DPNd 'AllSeg.mat']));
   SS(k,9)= num2cell(exist([DPNd 'AllSegCut.mat']));
    
   
   %Create Show OUtput
   SOut(k,1)=cellstr(Kdir(k).name);
   if exist([TPN 'Cell.mat'])
       load([TPN 'Cell.mat'])
       SOut(k,2)=cellstr(Cell.Age);
       SOut(k,3)=cellstr(Cell.Type);   
   end %end if cell exists 
       
   % Record Results
   if exist([DPNd 'Results.mat'])
       load([DPNd 'Results.mat'])
       SOut(k,4)=num2cell(Results.CellStats.AverageDensity);       
   end

end

SOut



