
clear all

%%Get directory with image folders
GPN=GetMyDir;
GPNd=dir(GPN); GPNd=GPNd(3:length(GPNd));
medSize = 3;  %  Define size of median filter kernal ( 1 for none)


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

        if medSize > 1
            'filtering',pause(.01)
            I=MyMedian(I,medSize);
        end
        Iname=Ifolder(1:find(Ifolder=='.',1)-1);
        if exist([GPN Iname])
            Iname = [Iname '_5'];
        end

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
if exist('Imaxs')
    Iwrite([GPN 'autoMaxTime'],Imaxs);
end
finished_files = ImageFolders






