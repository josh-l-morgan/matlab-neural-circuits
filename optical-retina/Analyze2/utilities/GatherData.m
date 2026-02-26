%Gather up info for smart guide to drop in a single folder

%KPN=uigetdir

%Drop=uigetdir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    TPN = [KPN '\' Kdir(k).name '\'];
    DPN = [TPN 'I\']   
     nTPN=[Drop '\' Kdir(k).name '\'];
     
    if ~exist([TPN 'find']) & exist([TPN 'Dots.mat'])
        mkdir([nTPN 'images'])
        copyfile([TPN 'Dots.mat'],[nTPN 'Dots.mat'])
        copyfile([TPN 'images\maxRaw.mat'],[nTPN 'images\maxRaw.mat'])       
        
    end
    pause(.1)
end