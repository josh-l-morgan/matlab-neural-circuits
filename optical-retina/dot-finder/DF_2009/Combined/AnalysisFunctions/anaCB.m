function[] = anaCB(TPN)

%Find cell body

colormap gray(255)

load([TPN 'ImageInfo.mat'])
xyum=ImageInfo{4};
zum=ImageInfo{5};


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'


load([TPN 'Dots.mat'])

%%{

TPNm = [TPN 'mask\'];
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));

if size(dTPNm,1)==1
    Ir=tiffread2([TPNm dTPNm(1).name]);
    IM=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
    for i = 1: size(Ir,2)
        IM(:,:,i)=Ir(i).data;
    end
    clear Ir
else

    clear I Ic
    IM(:,:)=imread([TPNm dTPNm(1).name]); %read
    IM(1,1,size(dTPNm,1))=0;
    c=0;
    for i=1:size(dTPNm,1)
        nam=dTPNm(i).name;
        naml=length(nam);
        if nam(naml-3:naml)=='.tif'
            c=c+1;
            IM(:,:,c)=imread([TPNm nam]);
        end
        PercentRead=i/size(dTPNm,1)*100
    end


end

[y x z] = ind2sub(size(IM),find(IM==3));
CBpos=[mean(y) mean(x) mean(z)];
Dots.Im.CBpos=CBpos;
%}

CBpos=Dots.Im.CBpos; %open and scale
CBpos(1:2)=CBpos(1:2)*xyum; CBpos(3)=CBpos(3)*zum;
Dpos=Dots.Pos; Dpos(:,1:2)=Dpos(:,1:2)*xyum; Dpos(:,3)=Dpos(:,3)*zum;

Dist2CB=dist(Dpos,CBpos);
Dots.Dist2CB=Dist2CB;
save([TPN 'Dots.mat'],'Dots')




