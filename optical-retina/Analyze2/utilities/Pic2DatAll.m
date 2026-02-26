

KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir),1)

%reads all the images in a folder, compiles them as a 3D matrix and saves
%the matrix as a mat file in the same directory that the folder is in

for k=1:size(Kdir,1)
    PPN=[KPN Kdir(k).name '\pics\']
    clear I 
    
    %Clear out that which is saved elsewhere
    if isdir([PPN 'BigCentroid']), rmdir([PPN 'BigCentroid'],'s'), end
    if isdir([PPN 'BigFilled']), rmdir([PPN 'BigFilled'],'s'), end
    if isdir([PPN 'BigIT']), rmdir([PPN 'BigIT'],'s'), end
    
    Pdir=dir(PPN); Pdir=Pdir(3:size(Pdir),1);
    for p = 1 : size(Pdir,1)
        Iname=Pdir(p).name;
        IPN=[PPN Iname];
        if isdir(IPN)
            Idir=dir(IPN); Idir=Idir(3:size(Idir,1));
            clear I
            I=imread([IPN '\' Idir(1).name]);
            [ys xs]=size(I);
            Clas=class(I);
            I=zeros(ys,xs,size(Idir,1),Clas);
            for i= 1:size(Idir,1)
                I(:,:,i)=imread([IPN '\' Idir(i).name]);
            end
            evalin('caller',[Iname '=' 'I' '; clear ' 'I']) %Switches file names
            save([PPN '\' Iname '.mat'],Iname);
            clear(Iname)
            rmdir(IPN,'s')
        end
    end
    
    
end
