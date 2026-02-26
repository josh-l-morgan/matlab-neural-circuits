%% Look at Data from all Cells while working out power issues

%% Read Data from Raid into Master List

clear all
%%Map Raid folders
Raid='\\128.208.64.36\wonglab\Josh\Analyzed\Output';
RaidDir=dir(Raid);
RD=struct2cell(RaidDir);
RD=RD(1,:); RD=RD(3:size(RD,2));

%%Compile Cell, Status, and Results into master list
for i=1:size(RD,2)
    if exist([Raid '\' char(RD(i)) '\Cell.mat']) 
    load([Raid '\' char(RD(i)) '\Cell.mat']) 
    Master(i).Cell=Cell; end
    if exist([Raid '\' char(RD(i)) '\Status.mat'])
     load([Raid '\' char(RD(i)) '\Status.mat']) 
    Master(i).Status=Status; end       
    if exist([Raid '\' char(RD(i)) '\data\Results.mat'])   
    load([Raid '\' char(RD(i)) '\data\Results.mat']) 
    Master(i).Results=Results; end        
end

%%pick out cells
for i=1:size(Master,2), c(i)=~isempty(Master(i).Cell); end
Master=Master(c); clear c

%%pick out results
for i=1:size(Master,2), c(i)=~isempty(Master(i).Results); end
Master=Master(c); clear c



%% Make lists of Cell types and Ages

ONs=logical(zeros(1,size(Master,2)));
OFFs=ONs; BIs=ONs; Others=ONs; Unknowns=ONs;

for i = 1:size(Master,2)
    
    %% read Cell type and mark appropriate matrix
    switch Master(i).Cell.Type
        case 'ON'
            ONs(i)=1;
        case 'OFF'
            OFFs(i)=1;
        case 'BI'
            BIs(i)=1;
        case 'Other'
            Others(i)=1;
        otherwise
            Unknowns(i)=1;
    end %End Cell type switch
    
    %% Put cell ages in Ages list
    Ages(i)=str2double(Master(i).Cell.Age);
    %% mark unknons as zero
    if isnan(Ages(i)), Ages(i)=0; end;
       
end

%%Extract Data
clear AverageDensity Arbor2DD Arbor21DD
for i = 1: size(Master,2)
   AverageDensity(i)=Master(i).Results.CellStats.AverageDensity; 
   Arbor1DD(i)=Master(i).Results.Arbor(1).DotDend;
   if size(Master(i).Results.Arbor,2)>1
       Arbor2DD(i) = Master(i).Results.Arbor(2).DotDend;
       MeanDD(i)=mean(Arbor1DD(i),Arbor2DD(i));
   else
       Arbor2DD(i)=0;
       MeanDD(i)=Arbor1DD(i);
   end
end

%% Display Results

hist(
        
        
        
        
        
        
        
        

