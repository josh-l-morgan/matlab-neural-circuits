
clear all

%%Get directory with image folders
GPN=GetMyDir;
GPNd=dir(GPN); GPNd=GPNd(3:length(GPNd));

%%Find image folders
ImageFolders={};
for i = 1:length(GPNd)
    Name=GPNd(i).name;
    siz=length(Name);
    if length(Name)>5;
    if Name(siz-5:siz)=='.files'
        ImageFolders(length(ImageFolders)+1,1)={Name};
    end
    end
end

%%write, max
clear Imaxs
TotalImages=length(ImageFolders)
for i = 1:length(ImageFolders)
   Report=['Running image ' num2str(i) ' of ' num2str(TotalImages)]
   
   'reading',pause(.01)
   Ifolder=cell2mat(ImageFolders(i));
   I=oifread([GPN Ifolder '\']);
   
   if size(size(I),2)>3
   
   'filtering',pause(.01)
   I=Shave(I,2); 
   Iname=Ifolder(1:find(Ifolder=='.',1)-1);
   
   'writing',pause(.01)
   I=juggleCh(I);
   Iwrite([GPN Iname],I)
   
   %% collect maxes
   Imax=max(I,[],4);
   if exist('Imaxs')
       [ys xs cs]=size(Imax);
       Imaxs(1:ys,1:xs,:,size(Imaxs,4)+1)=Imax;
   else
       Imaxs=Imax;
   end
  
   end
end

%%Write Max
Iwrite([GPN 'autoMaxTime'],Imaxs);

finished_files = ImageFolders






