%% Run all Programs up to 7/7/07

%{ 
Programs to run
Skeletonize from Mask
Dot Find



%}


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
        if ~isfield(SKStatus,'Ra'),SKStatus.Ra=0,end
        if ~isfield(SKStatus,'Ma'),SKStatus.Ma=0,end
        if ~isfield(SKStatus,'Rd'),SKStatus.Rd=0,end
        if ~isfield(SKStatus,'SG'),SKStatus.SG=0,end
        
        SKBusy=clock; save([TPN 'SKBusy.mat'],'SKBusy');
        SKBegin=clock; save([TPN 'SKBegin.mat'],'SKBegin')
        
        SKStatus
        %% Run Dot Processing
        if ~SKStatus.Sk 
            'Skeletonization'
            anaSkM(TPN,DPN)
            SKStatus.Sk=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.DF 
            'Finding Dots'
            dotFinder(TPN, DPN)
            SKStatus.DF=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Ra
            'Ratioing'
            anaRaM(TPN, DPN)
            SKStatus.Ra=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Ma
            'Masking'
            anaMaM(TPN, DPN)
            anaCB(TPN, DPN)
            anaFS(TPN, DPN) %check for shifts
            SKStatus.Ma=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end
        if ~SKStatus.Rd
            'Rounding'
            anaRd(TPN, DPN)
            SKStatus.Rd=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end
        
        if 1%~SKStatus.SG
            'Running SG once'
            anaSG(TPN, DPN)
 
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end

        delete([TPN 'SKBusy.mat'])
        SKEnd=1; save([TPN 'SKEnd.mat'],'SKEnd')
        clear SKStatus
    end %if necessary files are available
    

