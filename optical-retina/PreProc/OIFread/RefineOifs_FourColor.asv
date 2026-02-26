
clear all

%%Get directory with image folders
subFolders = recurFolder;

%% Get root
rootD = subFolders{1};
%mkdir([rootD '\imageIndex'])

medSize = 2;  %  Define size of median filter kernal ( 1 for none)

  
%%Find image folders
ImageFolders={};
for i = 1:length(subFolders)
    Name=subFolders{i};
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
    Report=['Running ' ImageFolders{i} '. Image ' num2str(i) ' of ' num2str(TotalImages)]
%         if ~exist([GPN Iname])

    'reading',pause(.01)
    Ifolder=cell2mat(ImageFolders(i));
    I=oifread([Ifolder '\']);
    
    if size(size(I),2)>3

        if medSize > 1
            'filtering',pause(.01)
            I=MyMedian(I,medSize);
        end


        %
%         %%Scale to Max
%         for i = 1:size(I,3)
%             minc=min(min(min(I(:,:,i,:))));
%             I(:,:,i,:)=I(:,:,i,:)-minc;
%             maxc=max(max(max(I(:,:,i,:))));
%             I(:,:,i,:)=double(I(:,:,i,:))*255/maxc;
%         end
        %%Scale to 8bit
        I = I/16;                   
        I=uint8(I);




        Iname=Ifolder(1:find(Ifolder=='.',1)-1);
        'writing',pause(.01)
        Ij=juggleCh(I,[1 2; 3 3; 4 1]);
        Iwrite(Iname,Ij)
        Ij=juggleCh(I,[1 2; 2 3; 4 1]);
        Iwrite([Iname '_both'],Ij)
        
%         end



    end

%     %%      collect maxes
%         Imax=max(I,[],4);
%         ID = Iname(length(rootD) + 2:length(Iname));
%         slashes= find(ID == '\');
%         if ~isempty(slashes)
%             ID(slashes) = 'X';
%         end
%         imwrite(Imax*2, [rootD '\imageIndex\' ID '.tif'],'Compression','none');
%     

    
%%Write Max

end
finished_files = ImageFolders

clear 




