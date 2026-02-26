function param=pipeline_getneighbortransabsrot(param)
%A function to compute the translations between neighboring tiles in a multi-tile stack
%This assumes that the rotation between tiles is negligible
%For use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 8th 2009

lparam=param.abstransneighbors;
lowfilt=lparam.lowestpixpercyc;
highfilt=lparam.highestpixpercyc; %5;

if lparam.docomputation
  
  transx_right=zeros(param.nrofslices, param.nrofrows, param.nrofcolumns-1);
  transy_right=zeros(param.nrofslices, param.nrofrows, param.nrofcolumns-1);
  transcorr_right=zeros(param.nrofslices, param.nrofrows, param.nrofcolumns-1);
  transx_down=zeros(param.nrofslices, param.nrofrows-1, param.nrofcolumns);
  transy_down=zeros(param.nrofslices, param.nrofrows-1, param.nrofcolumns);
  transcorr_down=zeros(param.nrofslices, param.nrofrows-1, param.nrofcolumns);
  
  for slice=1:1:param.nrofslices
    %compute translation vector to tile at the right
    if (param.nrofcolumns>1)
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns-1
          %load tile at (slice,row,column) as img1
          %here, the number of tiles is necessarily > 1
          name=sprintf(param.basescaledname,slice,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img1=imread(filename);
          img1=double(img1)/255;
          
          %load tile at (slice,row,column+1) as img2
          name=sprintf(param.basescaledname,slice,row,column+1);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img2=imread(filename);
          img2=double(img2)/255;
          
          [dx,dy,dc]=findtranslation(img1,img2,lowfilt,highfilt,30);
          transx_right(slice,row,column)=dx;
          transy_right(slice,row,column)=dy;
          transcorr_right(slice,row,column)=dc;
        end;
      end;
    end;
    
    %compute translation vector to tile below
    if (param.nrofrows>1)
      for row=1:1:param.nrofrows-1
        for column=1:1:param.nrofcolumns
          %load tile at (slice,row,column) as img1
          %here, the number of tiles is necessarily > 1
          name=sprintf(param.basescaledname,slice,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img1=imread(filename);
          img1=double(img1)/255;
          
          %load tile at (slice,row+1,column) as img2
          name=sprintf(param.basescaledname,slice,row+1,column);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img2=imread(filename);
          img2=double(img2)/255;
          
          [dx,dy,dc]=findtranslation(img1,img2,lowfilt,highfilt,30);
          transx_down(slice,row,column)=dx;
          transy_down(slice,row,column)=dy;
          transcorr_down(slice,row,column)=dc;
        end;
      end;
    end;
  end;
  
  param.abstransneighbors.transx_right=transx_right;
  param.abstransneighbors.transy_right=transy_right;
  param.abstransneighbors.transcorr_right=transcorr_right;
  param.abstransneighbors.transx_down=transx_down;
  param.abstransneighbors.transy_down=transy_down;
  param.abstransneighbors.transcorr_down=transcorr_down;
else
  %else just load in the estimates...
  transx_right=param.abstransneighbors.transx_right;
  transy_right=param.abstransneighbors.transy_right;
  transcorr_right=param.abstransneighbors.transcorr_right;
  transx_down=param.abstransneighbors.transx_down;
  transy_down=param.abstransneighbors.transy_down;
  transcorr_down=param.abstransneighbors.transcorr_down;  
end;

figure(30);
count=1;
if (param.nrofcolumns>1)
  for row=1:1:param.nrofrows
    for column=1:1:param.nrofcolumns-1
      subplot(param.nrofrows,param.nrofcolumns-1,count);
      plot(squeeze(transx_right(:,row,column)));
      hold on;
      plot(squeeze(transy_right(:,row,column)),'r');
      hold off;
      grid on;
      txt=sprintf('X (blue) and Y (red) for RIGHT vector from R%dC%d',row,column);
      title(txt);
      xlabel('Slice No.');
      ylabel('Coordinate value');
      count=count+1;
    end;
  end;
end;

figure(31);
count=1;
if (param.nrofrows>1)
  for row=1:1:param.nrofrows-1
    for column=1:1:param.nrofcolumns
      subplot(param.nrofrows-1,param.nrofcolumns,count);
      plot(squeeze(transx_down(:,row,column)));
      hold on;
      plot(squeeze(transy_down(:,row,column)),'r');
      hold off;
      grid on;
      txt=sprintf('X (blue) and Y (red) for DOWN vector from R%dC%d',row,column);
      title(txt);
      xlabel('Slice No.');
      ylabel('Coordinate value');
      count=count+1;
    end;
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test rigid alignment consistency

if (param.nrofrows>1)&&(param.nrofcolumns>1)
  %Alignment consistency can only be checked if there is a 2D-arrangement
  %of tiles
  
  transcon=zeros(param.nrofslices,2);
  transv=zeros(param.nrofslices,1);
  
  for slice=1:1:param.nrofslices
    dx1=transx_right(slice,1,1)+transx_down(slice,1,2); %dxy_r1c1_r1c2(wafer,slice,1)+dxy_r1c2_r2c2(wafer,slice,1);
    dy1=transy_right(slice,1,1)+transy_down(slice,1,2); %dxy_r1c1_r1c2(wafer,slice,2)+dxy_r1c2_r2c2(wafer,slice,2);
    dx2=transx_down(slice,1,1)+transx_right(slice,2,1); %dxy_r1c1_r2c1(wafer,slice,1)+dxy_r2c1_r2c2(wafer,slice,1);
    dy2=transy_down(slice,1,1)+transy_right(slice,2,1); %dxy_r1c1_r2c1(wafer,slice,2)+dxy_r2c1_r2c2(wafer,slice,2);
    transcon(slice,1)=dx2-dx1;
    transcon(slice,2)=dy2-dy1;
    transv(slice)=sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1));
  end;
  
  
  figure(32);
  %xv=squeeze(transcon(:,1));
  %yv=squeeze(transcon(:,2));
  %v=sqrt(xv.*xv+yv.*yv); %compute euclidian length of difference vector
  plot(transv); %plot(v);
  hold on;
  plot([1 param.nrofslices],[lparam.consistencythreshold lparam.consistencythreshold],'r--');
  hold off;
  xlabel('Slice No.');
  ylabel('Alignment consistency');
  %axis([0 60 -600 600]);
  grid on;
  
  param.abstransneighbors.transcon=transcon;
  param.abstransneighbors.transv=transv;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Recompute alignment on higher-resolution images for those slices in which consistency failed
  
  refineslicelist=find(param.abstransneighbors.transv>param.abstransneighbors.consistencythreshold);
  
  for i=1:1:max(size(refineslicelist))
    slice=refineslicelist(i);
    disp(slice);
    %compute translation vector to tile at the right
    if (param.nrofcolumns>1)
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns-1
          %load tile at (slice,row,column) as img1
          %here, the number of tiles is necessarily > 1
          name=sprintf(param.baserawname,slice,row,column);
          filename=sprintf('%s%s',param.rawdir,name);
          disp(filename);
          img1=imread(filename);
          img1=double(img1)/255;
          
          %load tile at (slice,row,column+1) as img2
          name=sprintf(param.baserawname,slice,row,column+1);
          filename=sprintf('%s%s',param.rawdir,name);
          disp(filename);
          img2=imread(filename);
          img2=double(img2)/255;
          
          [dx,dy,dc]=findtranslation(img1,img2,lowfilt*param.downscale.scale,highfilt*param.downscale.scale,40);
          transx_right(slice,row,column)=dx/param.downscale.scale;
          transy_right(slice,row,column)=dy/param.downscale.scale;
          transcorr_right(slice,row,column)=dc;
        end;
      end;
    end;
    
    %compute translation vector to tile below
    if (param.nrofrows>1)
      for row=1:1:param.nrofrows-1
        for column=1:1:param.nrofcolumns
          %load tile at (slice,row,column) as img1
          %here, the number of tiles is necessarily > 1
          name=sprintf(param.baserawname,slice,row,column);
          filename=sprintf('%s%s',param.rawdir,name);
          disp(filename);
          img1=imread(filename);
          img1=double(img1)/255;
          
          %load tile at (slice,row+1,column) as img2
          name=sprintf(param.baserawname,slice,row+1,column);
          filename=sprintf('%s%s',param.rawdir,name);
          disp(filename);
          img2=imread(filename);
          img2=double(img2)/255;
          
          [dx,dy,dc]=findtranslation(img1,img2,lowfilt*param.downscale.scale,highfilt*param.downscale.scale,40);
          transx_down(slice,row,column)=dx/param.downscale.scale;
          transy_down(slice,row,column)=dy/param.downscale.scale;
          transcorr_down(slice,row,column)=dc;
        end;
      end;
    end;
  end;
  
end;

% lowfilt=0;
% highfilt=5;
% docomputation=0;
% 
% if docomputation==1
% 
% dxy_r1c1_r1c2=zeros(nrofwafers,maxnrofslices,3); %dx, dy, maxcorrval
% dxy_r1c1_r2c1=zeros(nrofwafers,maxnrofslices,3);
% dxy_r1c2_r2c2=zeros(nrofwafers,maxnrofslices,3);
% dxy_r2c1_r2c2=zeros(nrofwafers,maxnrofslices,3);
% 
% for wafer=1:1:nrofwafers
%   for slice=1:1:maxnrofslices
%     %Load all four images of one slice
%     row=1; column=1;
%     name=sprintf('..\\Allslices_raw_10x_Matlab\\hip29nm_W%.3dS%.3dR%dC%d_10x.png',wafer,slice,row,column);
%     if exist(name,'file')==2 %only if this slice exists
%       imr1c1=imread(name); %load slice image
%       imr1c1=double(imr1c1)/255;
%     end;
%     row=1; column=2;
%     name=sprintf('..\\Allslices_raw_10x_Matlab\\hip29nm_W%.3dS%.3dR%dC%d_10x.png',wafer,slice,row,column);
%     if exist(name,'file')==2 %only if this slice exists
%       imr1c2=imread(name); %load slice image
%       imr1c2=double(imr1c2)/255;
%     end;
%     row=2; column=1;
%     name=sprintf('..\\Allslices_raw_10x_Matlab\\hip29nm_W%.3dS%.3dR%dC%d_10x.png',wafer,slice,row,column);
%     if exist(name,'file')==2 %only if this slice exists
%       imr2c1=imread(name); %load slice image
%       imr2c1=double(imr2c1)/255;
%     end;
%     row=2; column=2;
%     name=sprintf('..\\Allslices_raw_10x_Matlab\\hip29nm_W%.3dS%.3dR%dC%d_10x.png',wafer,slice,row,column);
%     if exist(name,'file')==2 %only if this slice exists
%       imr2c2=imread(name); %load slice image
%       imr2c2=double(imr2c2)/255;
%     end;
%     
%     dx=0; dy=0;
%     if exist(name,'file')==2
%       %compute all translation vectors between neighboring slice images
%       [dx,dy,dc]=findtranslation(imr1c1,imr1c2,lowfilt,highfilt,20);
%       dxy_r1c1_r1c2(wafer,slice,1)=dx;
%       dxy_r1c1_r1c2(wafer,slice,2)=dy;
%       dxy_r1c1_r1c2(wafer,slice,3)=dc;
%       
%       [dx,dy,dc]=findtranslation(imr1c1,imr2c1,lowfilt,highfilt,20);
%       dxy_r1c1_r2c1(wafer,slice,1)=dx;
%       dxy_r1c1_r2c1(wafer,slice,2)=dy;    
%       dxy_r1c1_r2c1(wafer,slice,3)=dc;   
%       
%       [dx,dy,dc]=findtranslation(imr1c2,imr2c2,lowfilt,highfilt,20);
%       dxy_r1c2_r2c2(wafer,slice,1)=dx;
%       dxy_r1c2_r2c2(wafer,slice,2)=dy;
%       dxy_r1c2_r2c2(wafer,slice,3)=dc;
%       
%       [dx,dy,dc]=findtranslation(imr2c1,imr2c2,lowfilt,highfilt,20);
%       dxy_r2c1_r2c2(wafer,slice,1)=dx;
%       dxy_r2c1_r2c2(wafer,slice,2)=dy;     
%       dxy_r2c1_r2c2(wafer,slice,3)=dc; 
%     end;
%   end;
% end;
% 
% else
%   name=sprintf('withinslicetrans_%d_%d.mat',lowfilt,highfilt);
%   load(name);
% end;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(1);
% subplot(1,2,1);
% plot(squeeze(dxy_r1c1_r1c2(:,:,1))');
% grid on;
% title('X-translation R1C1->R1C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% axis([0 60 0 512]);
% subplot(1,2,2);
% plot(squeeze(dxy_r1c1_r1c2(:,:,2))');
% grid on;
% title('Y-translation R1C1->R1C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% %axis([0 60 0 512]);
% 
% 
% figure(2);
% subplot(1,2,1);
% plot(squeeze(dxy_r1c1_r2c1(:,:,1))');
% grid on;
% title('X-translation R1C1->R2C1 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% axis([0 60 0 512]);
% subplot(1,2,2);
% plot(squeeze(dxy_r1c1_r2c1(:,:,2))');
% grid on;
% title('Y-translation R1C1->R2C1 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% %axis([0 60 0 512]);
% 
% figure(3);
% subplot(1,2,1);
% plot(squeeze(dxy_r1c2_r2c2(:,:,1))');
% grid on;
% title('X-translation R1C2->R2C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% axis([0 60 0 512]);
% subplot(1,2,2);
% plot(squeeze(dxy_r1c2_r2c2(:,:,2))');
% grid on;
% title('Y-translation R1C2->R2C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% %axis([0 60 0 512]);
% 
% figure(4);
% subplot(1,2,1);
% plot(squeeze(dxy_r2c1_r2c2(:,:,1))');
% grid on;
% title('X-translation R2C1->R2C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% axis([0 60 0 512]);
% subplot(1,2,2);
% plot(squeeze(dxy_r2c1_r2c2(:,:,2))');
% grid on;
% title('Y-translation R2C1->R2C2 (all wafers)');
% xlabel('Slice No.');
% ylabel('X-translation');
% %axis([0 60 0 512]);
% 
% 
% 
% figure(11);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c1_r1c2(wafer,:,1)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('X-translation R1C1->R1C2');
%   grid on;
% end;
% 
% figure(12);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c1_r1c2(wafer,:,2)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('Y-translation R1C1->R1C2');
%   grid on;
% end;
% 
% figure(13);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c1_r2c1(wafer,:,1)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('X-translation R1C1->R2C1');
%   grid on;
% end;
% 
% figure(14);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c1_r2c1(wafer,:,2)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('Y-translation R1C1->R2C1');
%   grid on;
% end;
% 
% figure(15);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c2_r2c2(wafer,:,1)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('X-translation R1C2->R2C2');
%   grid on;
% end;
% 
% figure(16);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r1c2_r2c2(wafer,:,2)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('Y-translation R1C2->R2C2');
%   grid on;
% end;
% 
% figure(17);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r2c1_r2c2(wafer,:,1)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('X-translation R2C1->R2C2');
%   grid on;
% end;
% 
% figure(18);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   plot(squeeze(dxy_r2c1_r2c2(wafer,:,2)));
%   %axis([0 60 0 512]);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('Y-translation R2C1->R2C2');
%   grid on;
% end;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Test rigid alignment consistency
% 
% transcon=zeros(nrofwafers,maxnrofslices,2);
% transv=zeros(nrofwafers,maxnrofslices);
% 
% for wafer=1:1:nrofwafers
%   for slice=1:1:maxnrofslices
%     dx1=dxy_r1c1_r1c2(wafer,slice,1)+dxy_r1c2_r2c2(wafer,slice,1);
%     dy1=dxy_r1c1_r1c2(wafer,slice,2)+dxy_r1c2_r2c2(wafer,slice,2);
%     dx2=dxy_r1c1_r2c1(wafer,slice,1)+dxy_r2c1_r2c2(wafer,slice,1);
%     dy2=dxy_r1c1_r2c1(wafer,slice,2)+dxy_r2c1_r2c2(wafer,slice,2);
%     transcon(wafer,slice,1)=dx2-dx1;
%     transcon(wafer,slice,2)=dy2-dy1;
%     transv(wafer,slice)=sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1));
%   end;
% end;
% 
% figure(19);
% for wafer=1:1:nrofwafers
%   subplot(2,3,wafer);
%   xv=squeeze(transcon(wafer,:,1));
%   yv=squeeze(transcon(wafer,:,2));
%   v=sqrt(xv.*xv+yv.*yv); %compute euclidian length of difference vector
%   plot(v);
%   t=sprintf('Wafer %d',wafer);
%   title(t);
%   xlabel('Slice No.');
%   ylabel('Alignment consistency');
%   axis([0 60 -600 600]);
%   grid on;
% end;
% 
% if docomputation==1
%   name=sprintf('withinslicetrans_%d_%d.mat',lowfilt,highfilt);
%   save(name,'dxy_r1c1_r1c2','dxy_r1c1_r2c1','dxy_r1c2_r2c2','dxy_r2c1_r2c2');
% end;
% 
% mean(abs(transv(abs(transv)<10)))
% mean(abs(transv(abs(transv)>=10)))
% size(transv(abs(transv)>=10),1)
