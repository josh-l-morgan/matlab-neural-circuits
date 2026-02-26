function param=photostagemap_ovregionalign(param)
%This function optimizes the alignment of overview images for a specific target region
%This assumes that the images have been pre-aligned
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.ovregionalign;

% startangle=0; endangle=360; nrofangles=36;
% rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
% scale1=lparam.downscale1; %first rough alignment can be skipped
scale1=lparam.downscale1;
scale2=lparam.downscale2;

for w=param.processwafers
  stsec=param.alignoverviewimages.stsec{w};
  nrofslices=size(stsec,1); %nrofslices=size(param.generatestagemap.linearslicepos{w},1);
  gavail=zeros(nrofslices,nrofslices);
  gtransx=zeros(nrofslices,nrofslices);
  gtransy=zeros(nrofslices,nrofslices);
  grot=zeros(nrofslices,nrofslices);
  gcorr=zeros(nrofslices,nrofslices);
  gmatrix=zeros(nrofslices,nrofslices,3,3);

  %generate strip section list
%   stsec=zeros(nrofslices,2);
%   i=1;
%   nrofstrips=size(param.generatestagemap.stripslicepos{w},1);
%   for st=1:1:nrofstrips
%     nrofsections=size(param.generatestagemap.stripslicepos{w}{st},1);
%     for sec=1:1:nrofsections
%       stsec(i,1)=st; stsec(i,2)=sec;
%       i=i+1;
%     end;
%   end;
    
  for i=1:1:nrofslices-1
    %Read image 1
    imagename1=[param.overviewimages.targetdirectory sprintf(param.overviewimages.filenametemplate,w,stsec(i,1),stsec(i,2))];
    img1=double(imread(imagename1));
    %Transform by global alignment and crop region of interest
    %alignedimg=zeros(size(img1,1),size(img1,2));
    %minx=1; miny=1; maxx=size(img1,2); maxy=size(img1,1); twidth=size(img1,2); theight=size(img1,1);
    minx=lparam.regionx(1); miny=lparam.regiony(1); maxx=lparam.regionx(2); maxy=lparam.regiony(2); 
    twidth=(maxx-minx+1); theight=(maxy-miny+1); %for 1:1 scaling, no zoom
    rmat=param.ovabsolutealign.gabsmatrix{w}{i}; %gabsmatrix{i};
    cutimg1=renderaffineimage(img1,rmat,minx,miny,maxx,maxy,twidth,theight,'cubic');
    if lparam.dofiltering
      cutimg1=bandpass2(cutimg1,lparam.filtminpx,lparam.filtmaxpx);
    end;

    figure(lparam.startfigure);
    imagesc(cutimg1); colormap gray; axis square;

    ms1cutimg1=imresize(cutimg1,1/scale1)/255; %downscalemean(img1,scale)/255;
    ms1cutimg1=ms1cutimg1-mean(mean(ms1cutimg1));
    ms2cutimg1=imresize(cutimg1,1/scale2)/255; %downscalemean(img1,scale2)/255;
    ms2cutimg1=ms2cutimg1-mean(mean(ms2cutimg1));

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
      
      minx=lparam.regionx(1); miny=lparam.regiony(1); maxx=lparam.regionx(2); maxy=lparam.regiony(2); 
      twidth=(maxx-minx+1); theight=(maxy-miny+1); %for 1:1 scaling, no zoom
      rmat=param.ovabsolutealign.gabsmatrix{w}{j}; %gabsmatrix{i};
      cutimg2=renderaffineimage(img2,rmat,minx,miny,maxx,maxy,twidth,theight,'cubic');
      if lparam.dofiltering
        cutimg2=bandpass2(cutimg2,lparam.filtminpx,lparam.filtmaxpx);
      end;
      
      %figure(lparam.startfigure+1);
      %imagesc(image2); colormap gray; axis square;
      figure(lparam.startfigure+1);
      cimage=uint8(zeros(size(cutimg1,1),size(cutimg1,2),3));
      cimage(:,:,1)=cutimg1;
      cimage(:,:,2)=cutimg2;
      cimage(:,:,3)=cutimg2;
      imshow(cimage); axis square;
      
      ms1cutimg2=imresize(cutimg2,1/scale1)/255; %downscalemean(img2,scale)/255;
      ms1cutimg2=ms1cutimg2-mean(mean(ms1cutimg2));
      ms2cutimg2=imresize(cutimg2,1/scale2)/255; %downscalemean(img2,scale2)/255;
      ms2cutimg2=ms2cutimg2-mean(mean(ms2cutimg2));
       
%       %Compute alignment between img1 and img2
%       disp('Coarse ...');
%       tic;
%       [transx,transy,corrval]=getrottransarrays(msimg1,msimg2,rots,'linear',0);
%       toc
%       [C,D]=max(corrval);
%       bestrot=rots(D);
%       figure(lparam.startfigure+2);
%       plot(rots,corrval,'r');
      
      %refine
      spread=lparam.rotrange; %spread=(360/nrofangles);
      rots=-spread:0.5:spread;
      repeat=1;
      while repeat==1
        disp('  Refine 1 ...');
        tic;
         txt=sprintf('Source image: %i, target image: %i',i,j); disp(txt);
        [transx,transy,corrval]=getrottransarrays(ms1cutimg1,ms1cutimg2,rots,'cubic',0);
        toc
        [C,D]=max(corrval);
        bestrot=rots(D);
        if (D==1)||(D==size(rots,2))
          if D==1
            rots=rots-spread;
          else
            rots=rots+spread;
          end;
        else
          repeat=0;
        end;
      end;
      figure(lparam.startfigure+2);
      plot(rots,corrval,'g');
      grid on;
      
      %refine more
      spread=0.5;
      rots2=bestrot-spread:0.1:bestrot+spread;
      repeat=1;
      while repeat==1
        disp('  Refine 2...');
        tic;
         txt=sprintf('Source image: %i, target image: %i',i,j); disp(txt);
        [transx,transy,corrval]=getrottransarrays(ms2cutimg1,ms2cutimg2,rots2,'cubic',0);
        toc
        [C,D]=max(corrval);
        bestrot=rots2(D);
        if (D==1)||(D==size(rots2,2))
          if D==1
            rots2=rots2-spread; %shift refinement window left
          else
            rots2=rots2+spread; %shift refinement window right
          end;
        else
          repeat=0;
        end;
      end;
      figure(lparam.startfigure+3);
      %hold on;
      plot(rots2,corrval*(scale2*scale2)/(scale1*scale1),'b');
      %hold off;
      grid on;
      %input('Press return...');
      
      %Compute and store alignment parameters
      gtransx(i,j)=(transx(D)+1)*scale2;
      gtransy(i,j)=(transy(D)+1)*scale2;
      grot(i,j)=rots2(D); 
      gcorr(i,j)=C;
      rot=-grot(i,j)*pi/180; tx=gtransx(i,j); ty=gtransy(i,j);
      rigidmatrix=[[cos(rot) -sin(rot) tx]; [sin(rot) cos(rot) ty]; [0 0 1]];
      gtransx(j,i)=rigidmatrix(1,3);
      gtransy(j,i)=rigidmatrix(2,3); 
      grot(j,i)=-grot(i,j);
      gcorr(j,i)=gcorr(j,i);
      %gmatrix(i,j,:,:)=inv(rigidmatrix);
      %gmatrix(j,i,:,:)=rigidmatrix;
      gmatrix(i,j,:,:)=inv(rigidmatrix)*(param.ovabsolutealign.gabsmatrix{w}{j}*inv(param.ovabsolutealign.gabsmatrix{w}{i}));
      gmatrix(j,i,:,:)=rigidmatrix*(param.ovabsolutealign.gabsmatrix{w}{i}*inv(param.ovabsolutealign.gabsmatrix{w}{j}));
      
      %rigidmatrix=squeeze(param.alignoverviewimages.gmatrix{w}(j,i,:,:))*rigidmatrix;
      %rigidmatrix=param.ovabsolutealign.gabsmatrix{w}{j}*rigidmatrix;
      %rigidmatrix=inv(squeeze(param.alignoverviewimages.gmatrix{w}(i,j,:,:))*inv(rigidmatrix));
      
      if lparam.showaligned==1
        %Show alignment of local regions
        salignrgb=zeros(size(cutimg1,1),size(cutimg1,2),3);
        scimg1=(cutimg1-min(min(cutimg1)))/((max(max(cutimg1))-min(min(cutimg1))));
        minx=1; miny=1; maxx=size(cutimg1,2); maxy=size(cutimg1,1); twidth=size(cutimg1,2); theight=size(cutimg1,1);
        %image2=renderaffineimage(img2,param.ovabsolutealign.gabsmatrix{w}{i}*inv(rigidmatrix),minx,miny,maxx,maxy,twidth,theight,'linear');
        rimage2=renderaffineimage(cutimg2,inv(rigidmatrix),minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg2=(rimage2-min(min(rimage2)))/((max(max(rimage2))-min(min(rimage2))));
        salignrgb(:,:,1)=scimg1;
        salignrgb(:,:,2)=scimg2;
        salignrgb(:,:,3)=scimg2;
        figure(lparam.startfigure+4);
        imshow(salignrgb);
        txt=sprintf('Slice %i aligned to slice %i',j,i); title(txt);
        
        salignrgb2=zeros(size(cutimg1,1),size(cutimg1,2),3);
        scimg2=(cutimg2-min(min(cutimg2)))/((max(max(cutimg2))-min(min(cutimg2))));
        %image1=renderaffineimage(img1,inv(param.ovabsolutealign.gabsmatrix{w}{j})*rigidmatrix,minx,miny,maxx,maxy,twidth,theight,'linear');
        rimage1=renderaffineimage(cutimg1,rigidmatrix,minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg1=(rimage1-min(min(rimage1)))/((max(max(rimage1))-min(min(rimage1))));
        salignrgb2(:,:,1)=scimg1;
        salignrgb2(:,:,2)=scimg2;
        salignrgb2(:,:,3)=scimg2;
        figure(lparam.startfigure+5);
        imshow(salignrgb2);
        txt=sprintf('Slice %i aligned to slice %i',i,j); title(txt);
      end;
      
      if lparam.showaligned2==1  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Show alignment of whole images in relative alignment (second slice transformed into the coordinate system of first slice)
        alignrgb=zeros(size(img1,1),size(img1,2),3);
        scimg1=(img1-min(min(img1)))/((max(max(img1))-min(min(img1))));
        minx=1; miny=1; maxx=size(img1,2); maxy=size(img1,1); twidth=size(img1,2); theight=size(img1,1);
        %image2=renderaffineimage(img2,param.ovabsolutealign.gabsmatrix{w}{i}*inv(rigidmatrix),minx,miny,maxx,maxy,twidth,theight,'linear');
        rimg2=renderaffineimage(img2,inv(rigidmatrix)*(param.ovabsolutealign.gabsmatrix{w}{j}*inv(param.ovabsolutealign.gabsmatrix{w}{i})),minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg2=(rimg2-min(min(rimg2)))/((max(max(rimg2))-min(min(rimg2))));
        alignrgb(:,:,1)=scimg1;
        alignrgb(:,:,2)=scimg2;
        alignrgb(:,:,3)=scimg2;
        figure(lparam.startfigure+6);
        imshow(alignrgb);
        txt=sprintf('Slice %i aligned to slice %i',j,i); title(txt);
        
        alignrgb2=zeros(size(img1,1),size(img1,2),3);
        scimg2=(img2-min(min(img2)))/((max(max(img2))-min(min(img2))));
        %image1=renderaffineimage(img1,inv(param.ovabsolutealign.gabsmatrix{w}{j})*rigidmatrix,minx,miny,maxx,maxy,twidth,theight,'linear');
        rimg1=renderaffineimage(img1,rigidmatrix*(param.ovabsolutealign.gabsmatrix{w}{i}*inv(param.ovabsolutealign.gabsmatrix{w}{j})),minx,miny,maxx,maxy,twidth,theight,'linear');
        scimg1=(rimg1-min(min(rimg1)))/((max(max(rimg1))-min(min(rimg1))));
        alignrgb2(:,:,1)=scimg1;
        alignrgb2(:,:,2)=scimg2;
        alignrgb2(:,:,3)=scimg2;
        figure(lparam.startfigure+7);
        imshow(alignrgb2);
        txt=sprintf('Slice %i aligned to slice %i',i,j); title(txt);
      end;
      if lparam.stopeach==1
        input('Press Return...');
      end;
    end;
  end;
  
  param.ovregionalign.gavail{w}=gavail;
  param.ovregionalign.gtransx{w}=gtransx;
  param.ovregionalign.gtransy{w}=gtransy;
  param.ovregionalign.grot{w}=grot;
  param.ovregionalign.gcorr{w}=gcorr;
  param.ovregionalign.gmatrix{w}=gmatrix;
end;


