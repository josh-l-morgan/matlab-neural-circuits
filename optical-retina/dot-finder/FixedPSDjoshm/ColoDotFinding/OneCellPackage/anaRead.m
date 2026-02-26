function[] = anaRead(TPN)

colormap gray(255)

DPN=[TPN 'I\']
Idir=dir(DPN); Idir=Idir(3:length(Idir));
Inam={};
for i = 1:length(Idir);  %Get all tiffs
    nam=Idir(i).name;
    if nam(length(nam)-3:length(nam))=='.tif'
        Inam{length(Inam)+1}=nam;
    end
end

%% Get dims
Is=imread([DPN Inam{1}]);
Iclass=class(Is);
[ys xs cs]=size(Is);
zs=length(Inam);

%% read
Iraw=zeros(ys, xs, zs, cs,Iclass);
'reading'
for i = 1: zs
   nam=Inam{i}; 
   Iraw(:,:,i,:)=imread([DPN nam]);
end
'done reading'


Imax=squeeze(max(Iraw,[],3));
image(Imax)
for i = 1:size(Imax,3)    
    subplot(1,3,i)
    image(Imax(:,:,i)*(500/double(max(max(Imax(:,:,i))))))
end

%% Get image info

xyum=.0695;
zum=.3;



'getting image info', pause(.1)
prompt = {'dendrite channel: ', 'post synaptic channel', 'colocalizing channel:  ',...
    'xy resolution : ', 'z resolution :'};
title = 'Define image channels. (0 if not applicable)';
nLines = 1;

ImageInfo= inputdlg(prompt,title,nLines,{'3','1','2',num2str(xyum),num2str(zum)});

for i = 1: length(ImageInfo)
    ImageInfo{i}=str2num(ImageInfo{i});
end

DenCh=ImageInfo{1};
PostCh=ImageInfo{2};
ColoCh=ImageInfo{3};
xyum=ImageInfo{4};
zum=ImageInfo{5};

save([TPN 'ImageInfo.mat'],'ImageInfo')

%% Write Channels


if DenCh,  Dend=Iraw(:,:,:,DenCh); save([TPN 'Dend.mat'],'Dend');end
if PostCh, Post=Iraw(:,:,:,PostCh); save([TPN 'Post.mat'],'Post'); end
if ColoCh, Colo=Iraw(:,:,:,ColoCh); save([TPN 'Colo.mat'],'Colo'); end















