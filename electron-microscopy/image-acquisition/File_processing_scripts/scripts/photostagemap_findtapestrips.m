function param=photostagemap_findtapestrips(param)
%This function is used to identify tape strips on a photographic image of a wafer.
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.tapestrips;

for w=param.processwafers
  waferradius=param.waferradius(w);
  stripdetectionthreshold=lparam.stripdetectionthreshold(w);
  %Load normalized image
  waferimagename=[param.normalizeddir sprintf(param.normalize.wafercroppedimagename,w)];
  txt=sprintf('Loading normalized wafer image %s ...',waferimagename); disp(txt);
  ccregistered=imread(waferimagename);
  
  figure(lparam.startfigure);
  imshow(ccregistered); axis equal;


  %Compute cropped wafer mask
  cwafermask=zeros(size(ccregistered,1),size(ccregistered,2));
  sqrad=(waferradius-lparam.wafermaskreduce)*(waferradius-lparam.wafermaskreduce);
  for y=1:1:size(cwafermask,1)
    dy=waferradius-y+1;
    for x=1:1:size(cwafermask,2)
      dx=waferradius-x+1;
      sqdist=dx*dx+dy*dy; %compute square of distance to center
      if sqdist<sqrad
        cwafermask(y,x)=1;
      end;
    end;
  end;
  cwafermask(end-lparam.wafermasklower:end,:)=0; %mask out lower end of image as well

  %Compute average brightness of wafer area
  cwmean=0; cwcount=0;
  for y=1:1:size(cwafermask,1)
    for x=1:1:size(cwafermask,2)
      if cwafermask(y,x)==1
        cwmean=cwmean+double(ccregistered(y,x));
        cwcount=cwcount+1;
      end;
    end;
  end;
  cwmean=cwmean/cwcount;

  %Mask out background
  cbwregistered=squeeze(mean(ccregistered,3)).*cwafermask;
  cbwregistered=cbwregistered+(1-cwafermask)*cwmean;
  figure(lparam.startfigure+1);
  imagesc(cbwregistered); colormap gray; axis equal;
  
  %Compute gradient image
  cbwdi=imageabsderivative(cbwregistered);
  figure(lparam.startfigure+2);
  imagesc(cbwdi); colormap gray; axis equal;


  %%%%%%%%%%%%%%%%%%%%% Compute Fourier Transform to estimate approximate orientation of strips
  disp('Computing Fourier transform and radializing..');
  lfreq=0.1; hfreq=0.5;
  
  %make divisible-by-two square image for fft
  if (size(cbwdi,1)~=size(cbwdi,2))||(mod(size(cbwdi,1),2)>0.5)
    dsize=min(size(cbwdi));
    if (mod(dsize,2)>0.5)%(Modulo value should be 0 or 1, thresholding it at 0.5 makes this rounding-error resistent)
      dsize=dsize-1;
    end;
  else
    dsize=size(im,1);
  end;
  cim=cbwdi(1:dsize,1:dsize); %crop image.
  fdim=abs(fftshift(fft2(cim))); %do fft
  figure(lparam.startfigure+3);
  imagesc(log(fdim));
  
  %do polar->cartesian transform using 0.1 degree resolution  
  hdsize=floor(dsize/2); %half of dsize
  lfreq=floor(lfreq*hdsize); %lowest frequency to consider
  hfreq=floor(hfreq*hdsize); %highest frequency to consider
  % sfdim=imresize(fdim,0.5);
  % isodd=floor(size(sfdim,1)/2)~=size(sfdim,1)/2;
  % if ~isodd
  %   sfdim=sfdim(2:end,2:end); dsize=dsize-1;
  % end;
  pdim=topolar4(fdim(2:dsize,2:dsize),([lfreq:hfreq])',([0:3599]*0.1*pi/180)','linear');
  %pdim=topolar4(sfdim(1:dsize/2,1:dsize/2),([lfreq:hfreq])',([0:3599]*0.1*pi/180)','linear');
  figure(lparam.startfigure+4);
  subplot(2,1,1);
  imagesc(log(pdim));
  
  pim2=sum(pdim); %integral over all relevant frequencies for each orientation
  pim2=pim2(1:1800)+pim2(1801:3600); %fold 0..360deg to 0..180deg)
  subplot(2,1,2);
  plot([0:1799]/10,pim2);
  grid on;
  [pwr,ang]=max(pim2);
  ang=-(ang-1)/10;
  
  %Compute angle in [-90,90]
  if (ang<-90)
    ang=ang+180;
  end;
  if (ang>90)
    ang=ang-180;
  end;
  txt=sprintf('Most prominent edge orientation in image: %.02f degrees.',ang); disp(txt);
  
  disp('Computing Radon transform ...');
  orient=[ang-lparam.stripangletolerance:lparam.radonanglestepping:ang+lparam.stripangletolerance];
  [rad,xp]=radon(cbwdi,orient);
  
  rad=rad((xp>=-waferradius)&(xp<=waferradius),:);
  xp=xp((xp>=-waferradius)&(xp<=waferradius));
  figure(lparam.startfigure+5);
  subplot(1,2,1);
  imagesc([orient(1) orient(end)],[1 size(rad,1)],rad);  
  xlabel('orientation (deg)');
  ylabel('image x-position');
  title('Result of Radon transform');
  
  if lparam.waferheightnormalize
    %Normalize edge probability by using wafer height
    lengthvec=real(sqrt(waferradius*waferradius-xp.*xp));
    lengthvec=lengthvec/max(lengthvec);
    figure(lparam.startfigure+6);
    subplot(1,2,1);
    plot(xp,lengthvec);
    grid on;
    
    norad=rad;
    norad(lengthvec>0,:)=norad(lengthvec>0,:)./((lengthvec(lengthvec>0))*ones(1,size(orient,2)));
    norad=norad*lparam.waferheightnormalizationfactor+rad*(1-lparam.waferheightnormalizationfactor);
    figure(lparam.startfigure+5);
    subplot(1,2,2);
    imagesc([orient(1) orient(end)],[1 size(norad,1)],norad);
    title('after wafer height normalization');
    
    figure(lparam.startfigure+6);
    subplot(1,2,2);
    plot(xp,mean(norad,2));
    grid on;
    rad=norad; %Write back to rad so that no distinction between rad and norad is necessary below
  end;
  
  if lparam.subtractminradon
    %Subtract the minimal value of the radon transform to improve the signal-to-noise ratio
    radmin=min(min(rad(lparam.edgemaskreduce:end-lparam.edgemaskreduce,:)));
    rad=rad-radmin; 
    rad=max(rad*0,rad);

    figure(lparam.startfigure+6);
    hold on;
    plot(xp,mean(rad,2),'r');
    hold off;
  end;

  %normalize
  radmax=max(max(rad));
  rad=rad/radmax;
  
  if lparam.radonthreshold>0
    %subtract additional threshold and normalize again
    rad=rad-lparam.radonthreshold;
    rad=max(rad*0,rad);
    radmax=max(max(rad));
    rad=rad/radmax;
  end;
  
  %Find parallel lines at 78 pixels distance
  rad2=sqrt(rad(lparam.tapewidth+1:end,:).*rad(1:end-lparam.tapewidth,:));
  mrad2=mean(rad2,2);
  %rad2=rad2.*rad2;
  figure(lparam.startfigure+7);
  subplot(1,2,1);
  imagesc([orient(1) orient(end)],[1 size(rad2,1)],rad2);
  subplot(1,2,2);
  plot(xp,mean(rad,2));
  hold on;
  plot(xp(1:size(mrad2,1)),mrad2,'g');
  hold off;
  grid on;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Find local maxima above stripdetectionthreshold
  nrlmax=0;
  for i=2:1:size(mrad2,1)-1
    if (mrad2(i)>stripdetectionthreshold)&&(mrad2(i)>mrad2(i-1))&&(mrad2(i)>mrad2(i+1))
      nrlmax=nrlmax+1;
    end;
  end;
  
  localmaxima=zeros(nrlmax,2);
  maxc=1;
  for i=2:1:size(mrad2,1)-1
    if (mrad2(i)>stripdetectionthreshold)&&(mrad2(i)>mrad2(i-1))&&(mrad2(i)>mrad2(i+1))
      localmaxima(maxc,1)=i;
      localmaxima(maxc,2)=mrad2(i);
      maxc=maxc+1;
    end;
  end
  
  figure(lparam.startfigure+7);
  subplot(1,2,2);
  hold on;
  for i=1:1:nrlmax
    plot(xp(localmaxima(i,1)),localmaxima(i,2),'r*');
  end;
  plot([xp(1) xp(end)],[stripdetectionthreshold stripdetectionthreshold],'r--');
  hold off;
  txt=sprintf('%d potential strips found.',nrlmax); disp(txt);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute location and orientation of kapton/peek strips
  localmaxima=[localmaxima [1:size(localmaxima,1)]' zeros(size(localmaxima,1))];
  if lparam.usegreedy==0
    % Count potential strips left-to-right
    %stripstart=zeros(nrlmax,1);
    localmaxima(1,4)=1; %stripstart(1)=1;
    if nrlmax>1
      for i=2:1:size(localmaxima,1)
        if (localmaxima(i)-localmaxima(i-1))>(0.5*lparam.tapewidth)  %assume that strips will not overlap by more than 50%.
          localmaxima(i,4)=1; %stripstart(i)=1;
        end;
      end;
    end;
  else
    % Find potential strips strong-to-weak
    %stripstart=zeros(nrlmax,1);
    localmaxima=sortrows(localmaxima,-2);
    for i=1:1:size(localmaxima,1)
      nearer=0;
      for j=1:1:size(localmaxima,1)
        d=localmaxima(i,1)-localmaxima(j,1);
        if (localmaxima(j,4)>0) && (abs(d)<(0.5*lparam.tapewidth)) %(stripstart(j)>0)
          nearer=1;
        end;
      end;
      if nearer==0
        localmaxima(i,4)=1; %stripstart(i)=1;
      end;
    end;
    localmaxima=sortrows(localmaxima,3);
  end;
  
  stripstart=localmaxima(:,4);
  
  txt=sprintf('%d probable strips found.',sum(stripstart)); disp(txt);
  pstripstart=localmaxima(find(stripstart),:);
  pstrip=zeros(size(pstripstart,1),2);
  pstrip(:,1)=xp(pstripstart(:,1)); %x-positions
  for i=1:1:size(pstripstart,1)
    [a,b]=max(rad2(pstripstart(i),:));
    pstrip(i,2)=orient(b); %edge orientation in degrees
  end;

  %Compute image regions of strips
  stripxcoords=zeros(size(pstrip,1),4);
  stripycoords=zeros(size(pstrip,1),4);
  radonctr=floor((size(cbwdi)+1)/2); %center of image as used by radon transform
  figure(lparam.startfigure);
  hold on;
  for i=1:1:size(pstrip,1)
    %   x=pstrip(i,1)*cos(pstrip(i,2)*pi/180)+radonctr(2);
    %   y=-pstrip(i,1)*sin(pstrip(i,2)*pi/180)+radonctr(1);
    length=sqrt(waferradius*waferradius-pstrip(i,1)*pstrip(i,1)); %cos(pstrip(i,1)/waferradius)*waferradius;
    k=length;
    x1=pstrip(i,1)*cos(pstrip(i,2)*pi/180)+k*sin(pstrip(i,2)*pi/180)+radonctr(2);
    y1=-pstrip(i,1)*sin(pstrip(i,2)*pi/180)+k*cos(pstrip(i,2)*pi/180)+radonctr(1);
    k=-length;
    x2=pstrip(i,1)*cos(pstrip(i,2)*pi/180)+k*sin(pstrip(i,2)*pi/180)+radonctr(2);
    y2=-pstrip(i,1)*sin(pstrip(i,2)*pi/180)+k*cos(pstrip(i,2)*pi/180)+radonctr(1);
    length=sqrt(waferradius*waferradius-(pstrip(i,1)+lparam.tapewidth)*(pstrip(i,1)+lparam.tapewidth)); %=cos((pstrip(i,1)+tapewidth)/waferradius)*waferradius;
    k=length;
    x3=(pstrip(i,1)+lparam.tapewidth)*cos(pstrip(i,2)*pi/180)+k*sin(pstrip(i,2)*pi/180)+radonctr(2);
    y3=(-pstrip(i,1)-lparam.tapewidth)*sin(pstrip(i,2)*pi/180)+k*cos(pstrip(i,2)*pi/180)+radonctr(1);
    k=-length;
    x4=(pstrip(i,1)+lparam.tapewidth)*cos(pstrip(i,2)*pi/180)+k*sin(pstrip(i,2)*pi/180)+radonctr(2);
    y4=(-pstrip(i,1)-lparam.tapewidth)*sin(pstrip(i,2)*pi/180)+k*cos(pstrip(i,2)*pi/180)+radonctr(1);
    plot([x1 x2 x4 x3 x1],[y1 y2 y4 y3 y1],'r');
    stripxcoords(i,:)=[x1 x2 x3 x4];
    stripycoords(i,:)=[y1 y2 y3 y4];
  end;
  hold off;
  
  param.tapestrips.nrofstrips(w)=size(pstrip,1);
  param.tapestrips.stripxcoords{w}=stripxcoords;
  param.tapestrips.stripycoords{w}=stripycoords;

  if (lparam.stopeach)&&(w<param.processwafers(end))
    input('Press Return...');
  end;
end;