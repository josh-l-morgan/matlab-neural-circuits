function param=photostagemap_normalize(param)
%This function contains routines to normalize a wafer photo (remove perspective).
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.normalize;

realredmarkerpos=[[-5.75 -5.75];[-5.75 5.75];[5.75 5.75];[5.75 -5.75]];
realgreenmarkerpos=[[-5.75 0];[0 5.75];[5.75 0];[0 -5.75]];

for w=param.processwafers
  waferradius=param.waferradius(w);
  imagename=[param.waferphotosdir sprintf(param.basewaferphotoname,w);];
  txt=sprintf('Reading image %s...',imagename); disp(txt);
  img=imread(imagename);

  disp('Looking for red markers...');
  redmap=zeros(size(img,1),size(img,2),3);
  redmap(:,:,1)=img(:,:,1)>lparam.redthres*img(:,:,2);
  redmap(:,:,2)=img(:,:,1)>lparam.redthres*img(:,:,3);
  redmap(:,:,3)=img(:,:,1)>30;
  rimg=redmap(:,:,1).*redmap(:,:,2).*redmap(:,:,3);
  
  [L,num] = bwlabeln(rimg);

  %keep only segments with more than 1000 pixels (relabel)
  txt=sprintf('%d red regions found. Selecting regions with more than 1000 pixels...',num); disp(txt);
  rnum=0;
  for i=1:1:num
    s=sum(sum(L==i));
    if (s>1000)
      rnum=rnum+1; L(L==i)=rnum;
    else
      rimg(L==i)=0; L(L==i)=0;
    end;
  end;

  if (rnum~=4)
    figure(lparam.startfigure);
    imagesc(rimg);
    txt=sprintf('ERROR: %d red marks found, instead of 4!',rnum); disp(txt);
    return;
  end;

  disp('Four red marks found.');

  %Compute coordinates of the red marks in the image
  redcenters=zeros(4,2);
  for i=1:1:4
    [x,y]=find(L==i);
    redcenters(i,1)=mean(x);
    redcenters(i,2)=mean(y);
  end;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  disp('Looking for green markers...');
  greenmap=zeros(size(img,1),size(img,2),3);
  greenmap(:,:,1)=img(:,:,2)>lparam.greenthres*img(:,:,1);
  greenmap(:,:,2)=img(:,:,2)>lparam.greenthres*img(:,:,3);
  greenmap(:,:,3)=img(:,:,2)>30;
  gimg=greenmap(:,:,1).*greenmap(:,:,2).*greenmap(:,:,3);

  [L,num] = bwlabeln(gimg);
  
  %keep only segments with more than 1000 pixels (relabel)
  txt=sprintf('%d green regions found. Selecting regions with more than 1000 pixels...',num); disp(txt);
  gnum=0;
  for i=1:1:num
    s=sum(sum(L==i));
    if (s>1000)
      gnum=gnum+1; L(L==i)=gnum;
    else
      gimg(L==i)=0; L(L==i)=0;
    end;
  end;
  
  if (gnum~=4)
    figure(lparam.startfigure);
    imagesc(gimg);
    txt=sprintf('ERROR: %d green marks found, instead of 4!',gnum); disp(txt);
    return;
  end;

  disp('Four green marks found.');

  %Compute coordinates of the green marks in the image
  greencenters=zeros(4,2);
  for i=1:1:4
    [x,y]=find(L==i);
    greencenters(i,1)=mean(x);
    greencenters(i,2)=mean(y);
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Sort red centers clockwise
  sredcenters=zeros(4,2);
  rcy=sortrows(redcenters,1); %sorted by Y coordinate
  rcx=sortrows(redcenters,2); %sorted by X coordinate
  %find upper-left and upper-right
  if ((rcy(1,1)==rcx(1,1))&&(rcy(1,2)==rcx(1,2))) || ((rcy(1,1)==rcx(2,1))&&(rcy(1,2)==rcx(2,2))) %if the first point is left
    sredcenters(1,:)=rcy(1,:); %upper left
    sredcenters(2,:)=rcy(2,:); %upper right
  else
    sredcenters(1,:)=rcy(2,:); %upper left
    sredcenters(2,:)=rcy(1,:); %upper right
  end;
  
  %find lower-right and lower-left
  if ((rcy(3,1)==rcx(3,1))&&(rcy(3,2)==rcx(3,2)))||((rcy(3,1)==rcx(4,1))&&(rcy(3,2)==rcx(4,2))) %if the third point is right
    sredcenters(3,:)=rcy(3,:); %lower right
    sredcenters(4,:)=rcy(4,:); %lower left
  else
    sredcenters(3,:)=rcy(4,:); %lower right
    sredcenters(4,:)=rcy(3,:); %lower left
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Sort green centers clockwise
  %Predictions of position of green markers: centers between red marker connecting lines
  p=[(sredcenters(1,:)+sredcenters(2,:))/2; (sredcenters(2,:)+sredcenters(3,:))/2; (sredcenters(3,:)+sredcenters(4,:))/2; (sredcenters(4,:)+sredcenters(1,:))/2;];
  sgreencenters=sortfit2d(greencenters,p); %Use distances to positions in p to sort greencenters
  figure(lparam.startfigure);
  imagesc(rimg*2+gimg);
  hold on;
  plot(sredcenters(:,2),sredcenters(:,1),'r*');
  plot(sredcenters(:,2),sredcenters(:,1),'r');
  plot(sgreencenters(:,2),sgreencenters(:,1),'k*');
  plot(sgreencenters(:,2),sgreencenters(:,1),'g');
  hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Compute coordinate transformation image -> paper
  % papermarkerpos=[realredmarkerpos*100; realgreenmarkerpos*100];
  % imagemarkerpos=[sredcenters; sgreencenters];
  papermarkerpos=[[realredmarkerpos(:,2) realredmarkerpos(:,1)]; [realgreenmarkerpos(:,2) realgreenmarkerpos(:,1)]]*10000/lparam.outputresolution;
  imagemarkerpos=[[sredcenters(:,2) sredcenters(:,1)]; [sgreencenters(:,2) sgreencenters(:,1)]];
  tform=cp2tform(imagemarkerpos,papermarkerpos,'projective');
  %tform=cp2tform(papermarkerpos,imagemarkerpos,'projective');
  [registered,xdata,ydata] = imtransform(img, tform);
  %itform = maketform('projective',tform.tdata.Tinv);
  %registered = imtransform(img, itform);
  figure(lparam.startfigure+1);
  imshow(registered);
  txt=sprintf('Registered photo at %.2f microns per pixel',lparam.outputresolution);
  title(txt);

  %Render cropped to marker centers
  mdist=5.75*10000/lparam.outputresolution;
  [cregistered,xdata,ydata] = imtransform(img, tform, 'XData',[-mdist mdist],'YData',[-mdist mdist]);
  figure(lparam.startfigure+2);
  imshow(cregistered);
  txt=sprintf('Registered cropped photo at %.2f microns per pixel',lparam.outputresolution);
  title(txt);
  
  %crop registered image to centers of color marks
  %compute location of upper left red marker in registered image
  outimagename=[param.normalizeddir sprintf(lparam.outimagename,w)];
  txt=sprintf('Writing cropped registered image to %s ...',outimagename); disp(txt);
  imwrite(cregistered,outimagename,lparam.outimagetype);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CROP TO WAFER
  
  %Find center of wafer
  bwregistered=squeeze(mean(cregistered,3));
  center=findwafer(bwregistered,waferradius,[575 575],50,lparam.startfigure+3); %Image is cropped to 575 pixels radius, search area is +-5mm
  txt=sprintf('Center of wafer %d is in the image at (x,y)=(%d,%d)',w,center(2),center(1)); disp(txt);
  param.normalize.wafercenter(w,:)=center;
  param.normalize.waferbboxy(w,:)=[center(1)-waferradius center(1)+waferradius];
  param.normalize.waferbboxx(w,:)=[center(2)-waferradius center(2)+waferradius];
  
  %Render wafer mask
  wafermask=zeros(size(bwregistered,1),size(bwregistered,2));
  sqrad=waferradius*waferradius;
  for y=1:1:size(wafermask,1)
    dy=center(1)-y;
    for x=1:1:size(wafermask,2)
      dx=center(2)-x;
      sqdist=dx*dx+dy*dy; %compute square of distance to center
      if sqdist<sqrad
        wafermask(y,x)=1;
      end;
    end;
  end;
  figure(lparam.startfigure+6);
  imagesc(wafermask); axis equal;
  wafermaskname=[param.normalizeddir sprintf(lparam.outmaskname,w)];
  txt=sprintf('Writing wafer mask image to %s ...',wafermaskname); disp(txt);
  imwrite(wafermask*255,wafermaskname,lparam.outimagetype);

  %Crop again to actual wafer
  ccregistered=cregistered([center(1)-waferradius:center(1)+waferradius],[center(2)-waferradius:center(2)+waferradius],:);
  figure(lparam.startfigure+7);
  imshow(ccregistered); axis equal;
  waferimagename=[param.normalizeddir sprintf(lparam.wafercroppedimagename,w)];
  title(waferimagename);
  txt=sprintf('Writing cropped registered image to %s ...',waferimagename); disp(txt);
  imwrite(ccregistered,waferimagename,lparam.outimagetype);
  
  if (lparam.stopeach)&&(w<param.processwafers(end))
    input('Press Return...');
  end;
end;