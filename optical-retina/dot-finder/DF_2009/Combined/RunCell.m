%% Run all Programs up to 1/19/09

%%Copyright Josh Morgan and Daniel Kerchensteiner 2005-2009


%{
Sets up sequence of programs to run to find puncta on labeled processes
%}


TPN = GetMyDir; %% Retrieves path for Folder that contains image folder (I)

%% Find out what has been done so far
if exist([TPN 'SKStatus.mat'])
    load([TPN 'SKStatus.mat']) %Retreive previous progress
else
    SKStatus.DF=1;  % Set up new status varible
end

%%Enter new Status fields if not already present
if ~isfield(SKStatus,'Read'),SKStatus.Read=1,end
if ~isfield(SKStatus,'Sk'),SKStatus.Sk=1,end
if ~isfield(SKStatus,'Ra'),SKStatus.Ra=1,end
if ~isfield(SKStatus,'Ma'),SKStatus.Ma=1,end
if ~isfield(SKStatus,'Rd'),SKStatus.Rd=1,end
if ~isfield(SKStatus,'SG'),SKStatus.SG=1,end


%% Find out what the user wants to do

%%Translate Status fields to user defaults
SKUser{1}=num2str(SKStatus.Read);
SKUser{2}=num2str(SKStatus.Sk);
SKUser{3}=num2str(SKStatus.DF);
SKUser{4}=num2str(SKStatus.Ra);
SKUser{5}=num2str(SKStatus.Ma);
SKUser{6}=num2str(1);

%%Set up input
'getting image info', pause(.1)
prompt = {'Read','Skeletonize ', 'Dot find',...
    'Ratio','Mask', 'Group'};
title = 'What analysis would you like to run';
nLines = 1;

%%Get user input
SKUser= inputdlg(prompt,title,nLines,{SKUser{1},SKUser{2},SKUser{3},...
    SKUser{4},SKUser{5},SKUser{6}});
pause(.1)

%%Translate input to status fields
SKStatus.Read= str2num(SKUser{1});
SKStatus.Sk = str2num(SKUser{2});
SKStatus.DF = str2num(SKUser{3});
SKStatus.Ra = str2num(SKUser{4});
SKStatus.Ma = str2num(SKUser{5});
SKStatus.Group = str2num(SKUser{6});


%% Set up necessary Directories
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps


%% Check out files
if SKStatus.Read
    anaRead(TPN)
    SKStatus.Read=0
    save([TPN 'SKStatus.mat'],'SKStatus')
end

%% Run Dot Processing
if SKStatus.Sk
    'Finding Skel'
    anaSk(TPN)
    SKStatus.Sk=0
    save([TPN 'SKStatus.mat'],'SKStatus')
end
if SKStatus.DF
    'Finding Dots'
    anaDF(TPN)
    SKStatus.DF=0
    save([TPN 'SKStatus.mat'],'SKStatus')
end
if SKStatus.Ra
    'Ratioing'
    anaRa(TPN)
    SKStatus.Ra=0
    save([TPN 'SKStatus.mat'],'SKStatus')
end
if SKStatus.Ma
    'Masking'
    anaMa(TPN)
    anaCB(TPN)
%     anaFSc(TPN, DPN) %check for shifts
    SKStatus.Ma=0
    save([TPN 'SKStatus.mat'],'SKStatus')
end

    'Running SG once'
    anaSG(TPN)
    
if SKStatus.Group 
    anaGroup(TPN)
    anaMakeUseOnec(TPN)
else
     anaMakeUseOne(TPN)
end
 


