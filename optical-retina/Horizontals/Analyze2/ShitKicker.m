

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
    M = exist([TPN 'data\D.mat']); %necessary for masking
    B = ~exist([TPN 'SKBegin.mat']); % has SK begun
    E = ~exist([TPN 'SKEnd.mat']); % Did SK finish
    
    if (I & R & M  & E) > 0 % if necessary files exist
        
        
        SKBegin=1; save([TPN 'SKBegin.mat'],'SKBegin')
        
        %% Run Dot Processing
        anaDF(TPN, DPN)
        anaRa(TPN, DPN)
        anaMa(TPN, DPN)
        anaRd
        
        SKEnd=1; save([TPN 'SKEnd.mat'],'SKEnd')

    end
    
    
end %Run all folders


