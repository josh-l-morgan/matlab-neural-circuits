function maskedI = axMask(TPN, I)
% TPN = GetMyDir;
TPNm = [TPN 'mask\'];
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));


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

%% Dilate mask
[ys xs zs] = size(D);
SE = strel('diamond', 3);
DD = zeros(ys,xs,zs);
for i = 1:zs
    DD(:,:,i) = imdilate(D(:,:,i),SE);
end
clear D


%% mask image
DD = uint8(DD);
DDrgb = uint8(zeros(size(I)));

for i = 1:3
    DDrgb(:,:,i,:) = DD;
end
clear DD
maskedI = I.*DDrgb;

clear I D

