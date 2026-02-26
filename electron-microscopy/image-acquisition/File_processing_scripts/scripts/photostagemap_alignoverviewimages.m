function param=photostagemap_alignoverviewimages(param)
%This function computes an alignment of a series of overview images taken
%at the locations predicted from the wafer photograph
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.alignoverviewimages;

%startangle=0; endangle=360; nrofangles=36;
startangle=lparam.level1startangle; endangle=lparam.level1endangle; nrofangles=lparam.level1nrofangles;
rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
scale1=lparam.downscale1;
scale2=lparam.downscale2;
scale3=lparam.downscale3;
%scale=64; scale2=16; scale3=8;

for w=param.processwafers
  nrofslices=sum(param.generatestagemap.linearslicepos{w}(:,3)==0); %size(param.generatestagemap.linearslicepos{w},1);
  
  gavail=zeros(nrofslices,nrofslices);
  gtransx=zeros(nrofslices,nrofslices);
  gtransy=zeros(nrofslices,nrofslices);
  grot=zeros(nrofslices,nrofslices);
  gcorr=zeros(nrofslices,nrofslices);
  gmatrix=zeros(nrofslices,nrofslices,3,3);

  %generate strip section list
  stsec=zeros(nrofslices,2);
  i=1;
  nrofstrips=size(param.generatestagemap.stripslicepos{w},1);
  for st=1:1:nrofstrips
    nrofsections=size(param.generatestagemap.stripslicepos{w}{st},1);
    for sec=1:1:nrofsections
      if (param.generatestagemap.stripslicepos{w}{st}(sec,3)==0) %only use 'good' sections for overview image alignment
        stsec(i,1)=st; stsec(i,2)=sec;
        i=i+1;
      end;
    end;
  end;
    
  for i=1:1:nrofslices-1
    imagename=[param.overviewimages.targetdirectory sprintf(param.overviewimages.filenametemplate,w,stsec(i,1),stsec(i,2))]
    %imagename=sprintf('./images/%02d.tif',i);
    img1=double(imread(imagename));
    figure(lparam.startfigure);
    imagesc(img1); colormap gray; axis square;
    msimg1=imresize(img1,1/scale1)/255; %downscalemean(img1,scale)/255;
    msimg1=msimg1-mean(mean(msimg1));
    ms2img1=imresize(img1,1/scale2)/255; %downscalemean(img1,scale2)/255;
    ms2img1=ms2img1-mean(mean(ms2img1));
    ms3img1=imresize(img1,1/scale3)/255; %downscalemean(img1,scale2)/255;
    ms3img1=ms3img1-mean(mean(ms3img1));
 
    %create target image array by using divisibility by powers of 2
    nr2p=1; v=i-1;
    while (floor(v/2)==v/2) && (i + 2^nr2p)<=nrofslices
      v=v/2; nr2p=nr2p+1;
    end;
    parr=zeros(nr2p,1);
    for p=1:1:nr2p
      parr(p)=i + 2^(p-1);
    end;
  
    for jc=1:1:size(parr,1)
      j=parr(jc);
      gavail(i,j)=1; gavail(j,i)=1;
      txt=sprintf('Source image: %i, target image: %i',i,j); disp(txt);
      %if i~=j %no need to align the slice with itself
      %imagename=sprintf('./images/%02d.tif',j);
      imagename2=[param.overviewimages.targetdirectory sprintf(param.overviewimages.filenametemplate,w,stsec(j,1),stsec(j,2))];
      img2=double(imread(imagename2));
      figure(lparam.startfigure+1);
      imagesc(img2); colormap gray; axis square;
      msimg2=imresize(img2,1/scale1)/255; %downscalemean(img2,scale)/255;
      msimg2=msimg2-mean(mean(msimg2));
      ms2img2=imresize(img2,1/scale2)/255; %downscalemean(img2,scale2)/255;
      ms2img2=ms2img2-mean(mean(ms2img2));
      ms3img2=imresize(img2,1/scale3)/255; %downscalemean(img2,scale2)/255;
      ms3img2=ms3img2-mean(mean(ms3img2));
      
      %Compute alignment between img1 and img2
      disp('Coarse ...');
      tic;
      [transx,transy,corrval]=getrottransarrays(msimg1,msimg2,rots,'linear',0);
      toc
      [C,D]=max(corrval);
      bestrot=rots(D);
      figure(lparam.startfigure+2);
      plot(rots,corrval,'r');
      
      %refine
      spread=((endangle-startangle)/nrofangles);
      rots2=bestrot-spread:0.5:bestrot+spread;
      repeat=1;
      while repeat==1
        disp(' Refine 1 ...');
        tic;
        [transx,transy,corrval]=getrottransarrays(ms2img1,ms2img2,rots2,'linear',0);
        txt=sprintf('Source image: %i, target image: %i',i,j); disp(txt);
        toc
        [C,D]=max(corrval);
        bestrot=rots2(D);
        if (D==1)||(D==size(rots2,2))
          if D==1
            rots2=rots2-spread;
          else
            rots2=rots2+spread;
          end;
        else
          repeat=0;
        end;
      end;
      figure(lparam.startfigure+2);
      hold on;
      plot(rots2,corrval*(scale2*scale2)/(scale1*scale1),'g');
      hold off;
      
      %refine more
      spread=0.5;
      rots3=bestrot-spread:0.1:bestrot+spread;
      repeat=1;
      while repeat==1
        disp('  Refine 2...');
        tic;
        [transx,transy,corrval]=getrottransarrays(ms3img1,ms3img2,rots3,'linear',0);
        txt=sprintf('Source image: %i, target image: %i',i,j); disp(txt);
        toc
        [C,D]=max(corrval);
        bestrot=rots3(D);
        if (D==1)||(D==size(rots3,2))
          if D==1
            rots3=rots3-spread; %shift refinement window left
          else
            rots3=rots3+spread; %shift refinement window right
          end;
        else
          repeat=0;
        end;
      end;
      figure(lparam.startfigure+3);
      %hold on;
      plot(rots3,corrval*(scale3*scale3)/(scale1*scale1),'b');
      %hold off;
      grid on;
      %input('Press return...');
      
      %Compute and store alignment parameters
      gtransx(i,j)=(transx(D)+1)*scale3;
      gtransy(i,j)=(transy(D)+1)*scale3;
      grot(i,j)=rots3(D); 
      gcorr(i,j)=C;
      rot=-grot(i,j)*pi/180; tx=gtransx(i,j); ty=gtransy(i,j);
      rigidmatrix=[[cos(rot) -sin(rot) tx]; [sin(rot) cos(rot) ty]; [0 0 1]];
      gtransx(j,i)=rigidmatrix(1,3);
      gtransy(j,i)=rigidmatrix(2,3); 
      grot(j,i)=-grot(i,j);
      gcorr(j,i)=gcorr(j,i);
      gmatrix(i,j,:,:)=inv(rigidmatrix);
      gmatrix(j,i,:,:)=rigidmatrix;
      
%       gtransx(i,j)=transx(D);
%       gtransy(i,j)=transy(D);
%       grot(i,j)=rots(D);
%       gcorr(i,j)=C;
      %end;
      
      if lparam.showaligned==1
        %Show alignment.
        alignrgb=zeros(size(img1,1),size(img1,2),3);
        scimg1=(img1-min(min(img1)))/((max(max(img1))-min(min(img1))));
        minx=1; miny=1; maxx=size(img1,2); maxy=size(img1,1); twidth=size(img1,2); theight=size(img1,1);
        image2=renderaffineimage(img2,inv(rigidmatrix),minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg2=(image2-min(min(image2)))/((max(max(image2))-min(min(image2))));
        alignrgb(:,:,1)=scimg1;
        alignrgb(:,:,2)=scimg2;
        alignrgb(:,:,3)=scimg2;
        figure(lparam.startfigure+4);
        imshow(alignrgb);
        txt=sprintf('Slice %i aligned to slice %i',j,i); title(txt);
        
        alignrgb2=zeros(size(img1,1),size(img1,2),3);
        scimg2=(img2-min(min(img2)))/((max(max(img2))-min(min(img2))));
        image1=renderaffineimage(img1,rigidmatrix,minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg1=(image1-min(min(image1)))/((max(max(image1))-min(min(image1))));
        alignrgb2(:,:,1)=scimg1;
        alignrgb2(:,:,2)=scimg2;
        alignrgb2(:,:,3)=scimg2;
        figure(lparam.startfigure+5);
        imshow(alignrgb2);
        txt=sprintf('Slice %i aligned to slice %i',i,j); title(txt);
      end;

    end;
    %figure(lparam.startfigure+4);
    % plot(squeeze(gcorr(i,:)));
    %grid on;
    %input('Press return...');
  end;
  
  % figure(6);
  % surf(gtransx);
  % figure(7);
  % surf(gtransy);
  % figure(8);
  % surf(grot);
  % figure(9);
  % surf(gcorr);
  
  %save('paramfile.mat','gavail','gtransx','gtransy','grot','gcorr','gmatrix');
  param.alignoverviewimages.gavail{w}=gavail;
  param.alignoverviewimages.gtransx{w}=gtransx;
  param.alignoverviewimages.gtransy{w}=gtransy;
  param.alignoverviewimages.grot{w}=grot;
  param.alignoverviewimages.gcorr{w}=gcorr;
  param.alignoverviewimages.gmatrix{w}=gmatrix;
  param.alignoverviewimages.stsec{w}=stsec;
end;


% gabsmatrix=cell(50);
% gabsmatrixcount=zeros(50,1);
% gabsmatrix{1}=[[1 0 0]; [0 1 0]; [0 0 1]]; gabsmatrixcount(1)=1;
% %Do recursive alignment to slice 1
% while sum(gabsmatrixcount==0)>0
% %for rec=1:1:nrofrecursions
%   for i=1:1:50
%     if gabsmatrixcount(i)>0 %source exists
%       for j=1:1:50
%         if gavail(i,j)
%           if gabsmatrixcount(j)==0
%             gabsmatrix{j}=gabsmatrix{i}*squeeze(gmatrix(i,j,:,:));
%             gabsmatrixcount(j)=1;
%           end;
%         end;
%       end;
%     end;
%   end;
% end;
% 
% %Render images
% for i=1:1:50
%   imagename=sprintf('./images/%02d.tif',i);
%   img1=double(imread(imagename));
%   alignedimg=zeros(size(img1,1),size(img1,2));
%   minx=1; miny=1; maxx=size(img1,2); maxy=size(img1,1); twidth=size(img1,2); theight=size(img1,1);
%   rmat=gabsmatrix{i};
%   image1=renderaffineimage(img1,rmat,minx,miny,maxx,maxy,twidth,theight,'linear');
%   alignedimg=(image1-min(min(image1)))/((max(max(image1))-min(min(image1))));
%   imagename=sprintf('./rigid/%02d.png',i);
%   imwrite(alignedimg,imagename,'png');
% end;