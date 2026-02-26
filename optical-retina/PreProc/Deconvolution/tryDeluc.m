%%PSF

TPN = GetMyDir
TPNd=dir(TPN); TPNd=TPNd(3:length(TPNd));
Name={};
for i = 1: length(TPNd)
    nam=TPNd(i).name;
    LN = length(nam);
    if nam(LN-3:LN) == '.tif'
        Name{length(Name)+1}=nam;
    end
end
for i = 1:length(Name)
    Rd(:,:,i,:)=imread([TPN '\' char(Name{i})]);
end
[ys xs zs cs] = size(I);
Rdmax=squeeze(max(Rd,[],3));
image(Rdmax),pause(.01)





TPN = GetMyDir
TPNd=dir(TPN); TPNd=TPNd(3:length(TPNd));
Name={};
for i = 1: length(TPNd)
    nam=TPNd(i).name;
    LN = length(nam);
    if nam(LN-3:LN) == '.tif'
        Name{length(Name)+1}=nam;
    end
end
for i = 1:length(Name)
    PSF(:,:,i,:)=imread([TPN '\' char(Name{i})]);
end
[ys xs zs cs] = size(I);
Pmax=squeeze(max(PSF,[],3));
image(Pmax),pause(.01)

DeLu=deconvlucy(Rd,PSF);




