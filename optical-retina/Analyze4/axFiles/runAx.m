%% Run all Programs up to 7/7/07

%{ 
Programs to run
Skeletonize from Mask
Dot Find



%}
clear; clc; clf;
xyum=.0412;
zum=.3;
itMin = 10;
dfofMin = 2;


TPN = GetMyDir;
DPN=[TPN 'I\']
    %%if files are available
    I = exist(DPN); %necessary for Dot Finder
    R = exist([TPN 'data\Threshold.mat']); %necessary for ratioing
    M = exist([TPN 'data\D.mat']) | exist([TPN 'mask']); %necessary for masking
    B = ~exist([TPN 'SKBegin.mat']); % has SK begun
    E = ~exist([TPN 'SKEnd.mat']); % Did SK finish
    Busy = ~exist([TPN 'SKBusy.mat']); % Is someone currently working on this cell
    
    if (I ) > 0 % if necessary files exist (or dont)
        
        if exist([TPN 'SKStatus.mat'])
            load([TPN 'SKStatus.mat'])
        else
            SKStatus.DF=0;
        end
        
        if ~isfield(SKStatus,'Sk'),SKStatus.Sk=0,end
        if ~isfield(SKStatus,'Co'),SKStatus.Co=0,end
        if ~isfield(SKStatus,'Ra'),SKStatus.Ra=0,end
        if ~isfield(SKStatus,'Ma'),SKStatus.Ma=0,end
        
        
        SKBusy=clock; save([TPN 'SKBusy.mat'],'SKBusy');
        SKBegin=clock; save([TPN 'SKBegin.mat'],'SKBegin')
        
        SKStatus
        %% Run Dot Processing
        if ~SKStatus.Sk 
            'Skeletonization'
            axSkel(TPN,xyum,zum);
            SKStatus.Sk=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Co
            'Finding Area'
            axConvex(TPN,xyum);
            SKStatus.Co=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        
        if ~SKStatus.DF 
            'Finding Dots'
            axDotNeuriteFinder(TPN,DPN)
            SKStatus.DF=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Ra
            'Ratioing'
            axAnaRaM(TPN, DPN)
            SKStatus.Ra=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if 1
            'Draw Dots'
            axDrawDots(TPN,DPN,xyum,zum)
            SKStatus.Ma=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end

        delete([TPN 'SKBusy.mat'])
        SKEnd=1; save([TPN 'SKEnd.mat'],'SKEnd')
        clear SKStatus
    end %if necessary files are available
    

