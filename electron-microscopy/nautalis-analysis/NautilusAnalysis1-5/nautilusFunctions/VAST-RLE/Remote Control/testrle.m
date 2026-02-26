%testrle.m

global vdata;

info=vdata.vast.getinfo();
miplevel=5;

% areaxmin=0;
% areaxmax=bitshift(info.datasizex,-miplevel)-1;
% areaymin=0;
% areaymax=bitshift(info.datasizey,-miplevel)-1;
% areazmin=0;
% areazmax=info.datasizez-1;
% 
% [segimageshell,res]=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,0,255,0,99,1);
% [segimagefull,res]=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,0,255,0,99,0);
% [segimageshell2,res]=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,16,255+16,16,99+16,1);
% [segimagefull2,res]=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,16,255+16,16,99+16,0);
% 
% figure(2);
% imagesc(squeeze(segimageshell(:,:,1))+squeeze(segimagefull(:,:,1)));
% figure(3);
% imagesc(squeeze(segimageshell(:,:,2))+squeeze(segimagefull(:,:,2)));
% figure(4);
% imagesc(squeeze(segimageshell(:,:,3))+squeeze(segimagefull(:,:,3)));
% 
% dshell=segimageshell2(1:end-16,1:end-16,1:end-16)-segimageshell(1:end-16,17:end,17:end);
% figure(5);
% imagesc(squeeze(dshell(:,:,2)));
% max(dshell(:))

% slice=3;
% 
% % img1f=vdata.vast.getsegimageRLEdecoded(miplevel,256,334,0,255,0,99,0);
% % img1=vdata.vast.getsegimageRLEdecoded(miplevel,256,334,0,255,0,99,1);
% % img2=vdata.vast.getsegimage(miplevel,256,334,0,255,0,99);
% % img2=reshape(img2,79,256,100);
% 
% img1f=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,256,401,0,99,0);
% img1=vdata.vast.getsegimageRLEdecoded(miplevel,0,255,256,401,0,99,1);
% img2=vdata.vast.getsegimage(miplevel,0,255,256,401,0,99);
% img2=reshape(img2,256,146,100);
% 
% figure(6);
% imagesc(squeeze(img1(:,:,1))+squeeze(img1f(:,:,slice)))
% figure(7);
% imagesc(squeeze(img1f(:,:,1))+squeeze(img2(:,:,slice)))
% figure(8);
% subplot(1,3,1);
% imagesc(squeeze(img1f(:,:,slice)));
% title('RLE full');
% subplot(1,3,2);
% imagesc(squeeze(img1(:,:,slice)));
% title('RLE hollow');
% subplot(1,3,3);
% imagesc(squeeze(img2(:,:,slice)));
% title('Raw full');

[segimage,values,numbers,bboxes,res] = vdata.vast.getsegimageRLEdecodedbboxes(5,0,255,0,255,0,100,0);

%compute bounding boxes from segimage to compare to bboxes
[v,n]=count_unique(segimage);

max(v-values)
max(double(numbers)-n)

bbx=zeros(max(size(v)),6);
for i=1:1:max(size(v))
  csx=squeeze(sum(sum(segimage==values(i),3),2));
  csy=squeeze(sum(sum(segimage==values(i),3),1));
  csz=squeeze(sum(sum(segimage==values(i),1),2));
  xmin=find(csx>0,1,'first'); xmax=find(csx>0,1,'last');
  ymin=find(csy>0,1,'first'); ymax=find(csy>0,1,'last');
  zmin=find(csz>0,1,'first'); zmax=find(csz>0,1,'last');
  bbx(i,:)=[xmin ymin zmin xmax ymax zmax];
end;