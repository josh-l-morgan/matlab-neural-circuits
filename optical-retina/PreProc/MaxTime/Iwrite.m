function[] = Iwrite(Path,I)

if exist(Path),rmdir(Path,'s'),end
   mkdir(Path) 

Slash=find(Path=='\');
name=Path(Slash(length(Slash))+1:length(Path));

Dims=size(size(I),2);
zs=size(I,Dims);
NumDig=fix(log10(zs))+1; % get digits
blank=zeros(NumDig,1);
blank=num2str(blank)';
% 
% %%Scale to 8bit
% for i = 1:3
%     minc=min(min(min(I(:,:,i,:))));
%     I(:,:,i,:)=I(:,:,i,:)-minc;
%     maxc=max(max(max(I(:,:,i,:))));
%     I(:,:,i,:)=double(I(:,:,i,:))*255/maxc;    
% end
% 
% I=uint8(I);

if Dims == 4
for i=1:zs
    Num=num2str(i);
    Numi=blank;
    Numi(size(Numi,2)-size(Num,2)+1:size(Numi,2))=Num;
    name1=[Path '\' name '_' Numi  '.tif'];
    imwrite(I(:,:,:,i),name1,'tif','Compression','none')
end

else
  for i=1:zs
    Num=num2str(i);
    Numi=blank;
    Numi(size(Numi,2)-size(Num,2)+1:size(Numi,2))=Num;
    name1=[Path '\' name '_' Numi  '.tif'];
    imwrite(I(:,:,i),name1,'tif','Compression','none')
end  
end