
KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir),1)

%reads all the images in a folder, compiles them as a 3D matrix and saves
%the matrix as a mat file in the same directory that the folder is in

for k=1:size(Kdir,1)
    PPN=[KPN Kdir(k).name '\pics\']
    if exist([PPN 'Irm.mat'])
        if ~isdir([PPN 'Irm'])

            if exist([PPN 'Irm.mat'])
                load([PPN 'Irm.mat']);
                imwriteNp([KPN Kdir(k).name '\'],Irm,'Irm')
                clear Irm
            end
        end
    end
end
