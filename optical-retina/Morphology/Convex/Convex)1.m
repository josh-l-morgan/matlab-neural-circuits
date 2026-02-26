clear all
TPN = GetMyDir;
TPNm = [TPN 'mask\'];
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));

if ~exist([TPN 'data']), mkdir([TPN 'data']),end
if ~exist([TPN 'temp']), mkdir([TPN 'temp']),end


%%Image Variables
xyum=.0695;
zum=.3;

prompt = {'XY voxel dimension: ','Z voxel dimension: '};
title = 'Image Info';
nLines = 1;

ImageInfo= inputdlg(prompt,title,nLines,{num2str(xyum),num2str(zum)});
xyum=str2num(ImageInfo{1});
zum=str2num(ImageInfo{2});

pause(.1)


aspect=zum/xyum;% ratio of z to xy dimentions
minObSize= 50; %% minimum size of object to be measured (default 20)
minFillSize = 10; %% Minimum size of continuous wave object (default 4)
maxSegLength = 5;  %% Maximum length of segments (default 2)

%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'reading image'


if size(dTPNm,1)==1
    Ir=tiffread2([TPNm dTPNm(1).name]);
    D=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
    for i = 1: size(Ir,2)
        D(:,:,i)=Ir(i).data;
    end
    clear Ir
else
    
    clear I Ic
    D(:,:)=imread([TPNm dTPNm(1).name]); %read
    D(1,1,size(dTPNm,1))=0;
    c=0;
    for i=1:size(dTPNm,1)
        nam=dTPNm(i).name;
        naml=length(nam);
        if nam(naml-3:naml)=='.tif'
            c=c+1;
            D(:,:,c)=imread([TPNm nam]);
        end
        PercentRead=i/size(dTPNm,1)*100
    end
    
   D=D==1;
 
end

%Find New sizes
[ys,xs,zs]=size(D);