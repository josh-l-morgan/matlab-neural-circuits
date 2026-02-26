function param=pipeline_readmetadatafromxml(param)
%Reads the image parameters from a XML metadata file as provided by the Lichtman lab.
%By Daniel Berger, July 13th 2009

param.getmetadata.sourceslice=param.reorder.sourceslice;
lparam=param.getmetadata;

param.getmetadata.tilewidth=zeros(1,param.nrofslices);
param.getmetadata.tileheight=zeros(1,param.nrofslices);
param.getmetadata.nrofrows=zeros(1,param.nrofslices);
param.getmetadata.nrofcolumns=zeros(1,param.nrofslices);
param.getmetadata.imgwidth_mu=zeros(1,param.nrofslices);
param.getmetadata.imgheight_mu=zeros(1,param.nrofslices);
param.getmetadata.pixelsize_nm=zeros(1,param.nrofslices);
param.getmetadata.xpos=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
param.getmetadata.ypos=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
param.getmetadata.rot=zeros(1,param.nrofslices);
%param.getmetadata.starttime=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);

for slice=1:param.nrofslices
  filename=sprintf(lparam.filenametemplate,lparam.sourceslice(slice));
  metadata=xml_readparamfromfile(filename);
  param.getmetadata.slicename{slice}=metadata.slicename;
  param.getmetadata.tilewidth(slice)=metadata.tilewidth;
  param.getmetadata.tileheight(slice)=metadata.tileheight;
  param.getmetadata.nrofrows(slice)=metadata.nrofrows;
  param.getmetadata.nrofcolumns(slice)=metadata.nrofcolumns;
  param.getmetadata.imgwidth_mu(slice)=metadata.imgwidth_mu;
  param.getmetadata.imgheight_mu(slice)=metadata.imgheight_mu;
  param.getmetadata.pixelsize_nm(slice)=metadata.tilewidth;
  param.getmetadata.xpos(slice,:,:)=metadata.xpos;
  param.getmetadata.ypos(slice,:,:)=metadata.ypos;
  param.getmetadata.rot(slice)=metadata.rot;
  param.getmetadata.starttime{slice,:,:}=metadata.starttime;
end;

figure(32);
plot(squeeze(param.getmetadata.xpos(:,1)),squeeze(param.getmetadata.ypos(:,1)));
hold on;
grid on;
plot(squeeze(param.getmetadata.xpos(:,1)),squeeze(param.getmetadata.ypos(:,1)),'o');
xlabel('xpos');
ylabel('ypos');
title('Stage position');
hold off;

figure(33);
plot(param.getmetadata.rot);
grid on;

rxpos=squeeze(param.getmetadata.xpos(:,1))'.*cos(param.getmetadata.rot*pi/180);
rxpos=rxpos-squeeze(param.getmetadata.ypos(:,1))'.*sin(param.getmetadata.rot*pi/180);
rypos=squeeze(param.getmetadata.xpos(:,1))'.*sin(param.getmetadata.rot*pi/180);
rypos=rypos+squeeze(param.getmetadata.ypos(:,1))'.*cos(param.getmetadata.rot*pi/180);
figure(34);
plot(rxpos,rypos);
hold on;
grid on;
plot(rxpos,rypos,'o');
xlabel('xpos');
ylabel('ypos');
title('Stage position');
hold off;

% pixsize=zeros(
% for slice=1:param.nrofslices
%   pixsize
% end;