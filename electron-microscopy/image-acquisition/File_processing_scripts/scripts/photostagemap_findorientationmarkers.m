function param=photostagemap_findorientationmarkers(param)
%This function identifies the location of the orientation markers 
%at the edge of the wafer in normalized images.
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.orientationmarkers;

for w=param.processwafers
  waferradius=param.waferradius(w);
  %Load normalized image
  waferimagename=[param.normalizeddir sprintf(param.normalize.wafercroppedimagename,w)];
  txt=sprintf('Loading normalized wafer image %s ...',waferimagename); disp(txt);
  ccregistered=imread(waferimagename);
  
  figure(lparam.startfigure);
  imshow(ccregistered); axis equal;
  
  %Find lower-left marker on wafer (next to yellow marker)
  %seek in region (24,648)-(78,708)
  winxmin=(-param.normalize.waferbboxx(w,1)+72+24); %*(waferradius/505); 
  winymin=(-param.normalize.waferbboxy(w,1)+90+648-40); %*(waferradius/505); %to adapt these to different wafer radiuses (for example 500)
  %yreg=mean(ccregistered(648:648+63,24:24+63,:),3);
  yreg=mean(ccregistered(winymin:winymin+63,winxmin:winxmin+63,:),3);
  figure(lparam.startfigure+1); %figure(9);
  colormap gray;
  subplot(2,3,1);
  imagesc(yreg);
  axis([1 64 1 64]); axis square;
  title('Original image (gray)');
  %clear lower left corner
  bgmask=zeros(size(yreg,1),size(yreg,2));
  for y=1:1:size(yreg,1)
    bgmask(y,1:floor(y/2.5)+4)=1;  %Add +2 here to go less towards the wafer edge
  end;  
  yrm=mean(mean(yreg(bgmask==0)));
%   for y=1:1:size(yreg,1)
%     yreg(y,1:floor(y/2.5)+4)=yrm;  %Add +2 here to go less towards the wafer edge
%   end;
  for y=1:1:size(yreg,1)
    for x=1:1:size(yreg,2)
      if bgmask(y,x)==1
        yreg(y,x)=yrm;
      end;
    end;
  end;
  subplot(2,3,2);
  imagesc(yreg);
  axis([1 64 1 64]); axis square;
  title('Wafer region');
  %Threshold at mid of brightness range to find marker
  thres=(max(max(yreg))+min(min(yreg)))/2;
  yregmask=yreg<thres;
  
  %Look for regions in thresholded image and find one which is large enough and closest to the expected location
  [L,num] = bwlabeln(yregmask);
  if num==1 %only one region found: take mean coordinate
    [i,j]=find(yregmask);
    ymarkerpos=[mean(i) mean(j)];
  else %more than one region found: choose region closest to expected coordinates
    regchoose=0; expx=32; expy=32; mind=1000;
    for i=1:1:num
      syr=(L==i);
      [yl,xl]=find(syr);
      if (size(yl,1)>50)  %only accept large enough regions
        myl=mean(yl);
        mxl=mean(xl);
        d=sqrt((expx-mxl)*(expx-mxl)+(expy-myl)*(expy-myl));
        if (d<mind)
          mind=d;
          regchoose=i;
        end;
      end;
    end;
    if (regchoose==0)
      disp('WARNING: Rotation marker 1 not found. Using background average');
    end;
    [yl,xl]=find(L==regchoose); %if all regions are too small, this will use all non-labeled pixels
    ymarkerpos=[mean(yl) mean(xl)];
  end;
  
  subplot(2,3,3);
  imagesc(yregmask);
  axis([1 64 1 64]); axis square;
  title('Marker 1 detected');
  
  gymarkerpos=[ymarkerpos(1)+winymin-1 ymarkerpos(2)+winxmin-1];
  figure(lparam.startfigure+1);
  subplot(2,3,1); hold on; plot(ymarkerpos(2),ymarkerpos(1),'r*'); hold off;
  subplot(2,3,2); hold on; plot(ymarkerpos(2),ymarkerpos(1),'r*'); hold off;
  subplot(2,3,3); hold on; plot(ymarkerpos(2),ymarkerpos(1),'r*'); hold off;
  figure(lparam.startfigure); hold on; plot(gymarkerpos(2),gymarkerpos(1),'r*'); hold off;
  
  %Find lower-left marker on wafer (next to blue marker)
  %seek in region
  %winxmin=880*(waferradius/505); winymin=200*(waferradius/505);
  winxmin=(-param.normalize.waferbboxx(w,1)+72+880); 
  winymin=(-param.normalize.waferbboxy(w,1)+90+200-30);
  breg=mean(ccregistered(winymin:winymin+63,winxmin:winxmin+63,:),3);
  figure(lparam.startfigure+1);
  colormap gray;
  subplot(2,3,4);
  imagesc(breg);
  axis([1 64 1 64]); axis square;
  title('Original image (gray)');
  %clear lower left corner
  bgmask=zeros(size(breg,1),size(breg,2));
  for y=1:1:size(breg,1)
    %bgmask(y,1:floor(y/2.5)+4)=1;  %Add +2 here to go less towards the wafer edge
    bgmask(y,floor(y/1.5)+20:end)=1;
  end;  
  brm=mean(mean(breg(bgmask==0)));
  for y=1:1:size(breg,1)
    for x=1:1:size(breg,2)
      if bgmask(y,x)==1
        breg(y,x)=brm;
      end;
    end;
  end;
  
%   brm=mean(mean(breg));
%   for y=1:1:size(breg,1)
%     breg(y,floor(y/1.5)+20:end)=brm;  %Add +2 here to go less towards the wafer edge
%   end;
  subplot(2,3,5);
  imagesc(breg);
  axis([1 64 1 64]); axis square;
  title('Wafer region');
  %Threshold at mid of brightness range to find marker
  thres=(max(max(breg))+min(min(breg)))/2;
  bregmask=breg<thres;
  
  %Look for regions in thresholded image and find one which is large enough and closest to the expected location
  [L,num] = bwlabeln(bregmask);
  if num==1 %only one region found: take mean coordinate
    [yl,xl]=find(bregmask);
    bmarkerpos=[mean(yl) mean(xl)];
  else %more than one region found: choose region closest to expected coordinates
    regchoose=0; expx=35; expy=35; mind=1000;
    for i=1:1:num
      sbr=(L==i);
      [yl,xl]=find(sbr);
      if (size(yl,1)>20)  %only accept large enough regions
        myl=mean(yl);
        mxl=mean(xl);
        d=sqrt((expx-mxl)*(expx-mxl)+(expy-myl)*(expy-myl));
        if (d<mind)
          mind=d;
          regchoose=i;
        end;
      end;
    end;
    if (regchoose==0)
      disp('WARNING: Rotation marker 2 not found. Using background average');
    end;
    [yl,xl]=find(L==regchoose); %if all regions are too small, this will use all non-labeled pixels
    bmarkerpos=[mean(yl) mean(xl)];
  end;
  
  subplot(2,3,6);
  imagesc(bregmask);
  axis([1 64 1 64]); axis square;
  title('Marker 2 detected');
  % [i,j]=find(bregmask); %now done above
  % bmarkerpos=[mean(i) mean(j)];
  gbmarkerpos=[bmarkerpos(1)+winymin-1 bmarkerpos(2)+winxmin-1];
  figure(lparam.startfigure+1);
  subplot(2,3,4); hold on; plot(bmarkerpos(2),bmarkerpos(1),'r*'); hold off;
  subplot(2,3,5); hold on; plot(bmarkerpos(2),bmarkerpos(1),'r*'); hold off;
  subplot(2,3,6); hold on; plot(bmarkerpos(2),bmarkerpos(1),'r*'); hold off;
  figure(lparam.startfigure); hold on; plot(gbmarkerpos(2),gbmarkerpos(1),'r*'); hold off;
  
  %Compute angle of ymarkerpos and bmarkerpos with respect to center of wafer.
  %The angle is in degrees, mathematically positive, starting from the positive x direction in the image
  ymdy=gymarkerpos(1)-waferradius; ymdx=gymarkerpos(2)-waferradius;
  bmdy=gbmarkerpos(1)-waferradius; bmdx=gbmarkerpos(2)-waferradius;
  yangle=(atan2(ymdy,-ymdx)+pi)*180/pi;
  bangle=(atan2(bmdy,-bmdx)+pi)*180/pi;
  txt=sprintf('Identified marker directions: Yellow: %.02f degrees, Blue: %.02f degrees',yangle, bangle); disp(txt);
  
  param.orientationmarkers.yangle(w)=yangle;
  param.orientationmarkers.bangle(w)=bangle;
  param.orientationmarkers.gbmarkerpos(w,1)=gbmarkerpos(1);
  param.orientationmarkers.gbmarkerpos(w,2)=gbmarkerpos(2);
  param.orientationmarkers.gymarkerpos(w,1)=gymarkerpos(1);
  param.orientationmarkers.gymarkerpos(w,2)=gymarkerpos(2);
  
  if (lparam.stopeach)&&(w<param.processwafers(end))
    input('Press Return...');
  end;
end;
