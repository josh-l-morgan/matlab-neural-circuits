

KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
if ~exist([KPN 'Combos']),mkdir([KPN 'Combos']),end

for i = 1:size(Kdir,1)
    TPN = [KPN '\' Kdir(i).name '\'];
    DPN = [TPN 'I\']        
    
    
    %%if files are available
    I = exist(DPN); %necessary for Dot Finder
    R = exist([TPN 'data\Threshold.mat']); %necessary for ratioing
    M = exist([TPN 'data\D.mat']) | exist([TPN 'mask']); %necessary for masking
    B = ~exist([TPN 'SKBegin.mat']); % has SK begun
    E = ~exist([TPN 'SKEnd.mat']); % Did SK finish
    Busy = ~exist([TPN 'SKBusy.mat']); % Is someone currently working on this cell
    
    if (I & R & M ) > 0 % if necessary files exist (or dont)
        
        if exist([TPN 'SKStatus.mat'])
            load([TPN 'SKStatus.mat'])
        else
            SKStatus.DF=0;
        end
        
        if ~isfield(SKStatus,'Ra'),SKStatus.Ra=0,end
        if ~isfield(SKStatus,'Ma'),SKStatus.Ma=0,end
        if ~isfield(SKStatus,'Rd'),SKStatus.Rd=0,end
        if ~isfield(SKStatus,'Crit'),SKStatus.Crit=0,end
        
        SKBusy=clock; save([TPN 'SKBusy.mat'],'SKBusy');
        SKBegin=clock; save([TPN 'SKBegin.mat'],'SKBegin')
        
        SKStatus
        %% Run Dot Processing
        if ~SKStatus.DF 
            'Finding Dots'
            anaDF(TPN, DPN)
            SKStatus.DF=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Ra
            'Ratioing'
            anaRa(TPN, DPN)
            SKStatus.Ra=1
            save([TPN 'SKStatus.mat'],'SKStatus')
        end
        if ~SKStatus.Ma
            'Masking'
            anaMa(TPN, DPN)
            SKStatus.Ma=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end
        if ~SKStatus.Rd
            'Rounding'
            anaRd(TPN, DPN)
            SKStatus.Rd=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end
        if ~SKStatus.Crit
            'Applying Criteria'
            anaCrit(TPN, DPN)
            SKStatus.Crit=1
            save([TPN 'SKStatus.mat'],'SKStatus')            
        end

        delete([TPN 'SKBusy.mat'])
        SKEnd=1; save([TPN 'SKEnd.mat'],'SKEnd')
        clear SKStatus
    end %if necessary files are available
    
    
end %Run all folders


