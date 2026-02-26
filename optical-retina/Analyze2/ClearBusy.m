%%Clear Busys
%%Remove busy signal from all cells



KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
if ~exist([KPN 'Combos']),mkdir([KPN 'Combos']),end

for i = 1:size(Kdir,1)
    TPN = [KPN '\' Kdir(i).name '\'];
    DPN = [TPN 'I\']        
    
    if exist([TPN 'SKBusy.mat'])
        delete([TPN 'SKBusy.mat'])
    end
    
end %Run all folders



