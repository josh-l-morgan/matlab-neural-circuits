function param=pipeline_refinerelrot(param)
%A function to compute relative rotations and translations between consecutive slices
%based on alignSMG.m
%For use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 3 2009

lparam=param.refinerelrot;

rotarr=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
transxarr=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
transyarr=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);

%lparam.corotarr contains the coarse rotation estimates to be used

for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    firstslice=1;
    for slice=2:1:param.nrofslices
      %Load image(s)
      if firstslice==1
        firstslice=0;
        if (param.nrofrows==1)&&(param.nrofcolumns==1)
          name=sprintf(param.basescaledname,slice-1);
        else
          name=sprintf(param.basescaledname,slice-1,row,column);
        end;
        filename=sprintf('%s%s',param.scaleddir,name);
        disp(filename);
        img1=imread(filename);
        %img1=bandpass2(img1,0,100);
        msimg1=downscalemean(img1,lparam.scale)/255;
        msimg1=msimg1-mean(mean(msimg1));
      else
        img1=img2;
        msimg1=msimg2;
      end;
      
      if (param.nrofrows==1)&&(param.nrofcolumns==1)
        name=sprintf(param.basescaledname,slice);
      else
        name=sprintf(param.basescaledname,slice,row,column);
      end;
      filename=sprintf('%s%s',param.scaleddir,name);
      disp(filename);
      img2=imread(filename);
      %img2=bandpass2(img2,0,100);
      msimg2=downscalemean(img2,lparam.scale)/255;
      msimg2=msimg2-mean(mean(msimg2));
      
      %startangle=0; endangle=360; nrofangles=lparam.nrofangles;
      %rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
      startangle=lparam.corotarr(slice,row,column)-lparam.range/2;
      endangle=lparam.corotarr(slice,row,column)+lparam.range/2;
      nrofangles=lparam.nrofangles;
      rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
      
      txt=sprintf('Computing maximal cross-correlation peak for all %d relative rotations ...',lparam.nrofangles);
      disp(txt);
      [transx,transy,corrval]=getrottransarrays(msimg1,msimg2,rots,'linear',0);
      [C,D]=max(corrval);
      rotarr(slice,row,column)=rots(D);
      if rotarr(slice,row,column)>180
        rotarr(slice,row,column)=rotarr(slice,row,column)-360;
      end;
      %corrarr(slice,:)=corrval;
      transxarr(slice,row,column)=(transx(D)+0.5)*lparam.scale;
      transyarr(slice,row,column)=(transy(D)+0.5)*lparam.scale;
      
      %visualize
      figure(1);
      imagesc(msimg1);
      colormap(gray);
      axis square;

      %s=floor(size(msimg2,1)*sqrt(2)/2);
      %rmsimg2=rotateimage(msimg2,rots(D)*pi/180,s,s,'linear');
      %         figure(2);
      %         imagesc(rmsimg2);
      %         colormap(gray);
      %         axis square;
      
      s=floor(size(img2,1)*sqrt(2)/2);
      figure(3);
      x1=(size(img1,2)/2)-(s/2)+transxarr(slice); %transx is horizontal, transy is vertical.
      y1=(size(img1,1)/2)-(s/2)+transyarr(slice);
      x2=min(size(img1,2),x1+s);
      y2=min(size(img1,1),y1+s);
      x1=max(x1,1); y1=max(y1,1); %lower image clipping
      imagesc(img1(y1:y2,x1:x2));
      axis([1 s 1 s]);
      colormap(gray);
      axis square;
      txt=sprintf('Slice %d, which is compared with Slice %d',slice-1,slice);
      title(txt);
      
      figure(4);
      rimg2=rotateimage(double(img2),rots(D)*pi/180,s,s,'linear');
      imagesc(rimg2);
      colormap(gray);
      axis square;
      txt=sprintf('Slice %d, which is compared with Slice %d',slice,slice-1);
      title(txt);
      
      pause(1);
      if lparam.stopeach
        input('Return..');
      end;
    end;
    
    figure(10+(row-1)*param.nrofcolumns+(column-1));
    subplot(3,1,1);
    plot(squeeze(rotarr(:,row,column)),'r');
    grid on;
    title('Relative Rotation');
    xlabel('Slice');
    
    subplot(3,1,2);
    plot(squeeze(transxarr(:,row,column)),'r');
    grid on;
    title('Relative X Translation');
    xlabel('Slice');
    
    subplot(3,1,3);
    plot(squeeze(transyarr(:,row,column)),'r');
    grid on;
    title('Relative Y Translation');
    xlabel('Slice');
    
  end;
end;

param.refinerelrot.rotarr2=rotarr;
param.refinerelrot.transxarr2=transxarr;
param.refinerelrot.transyarr2=transyarr;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fine-tune orientation estimate
% 
% compute3=0;
% scale=4;
% 
% if compute3==1
%   rotarr2=zeros(nrofslices,nrofrows,nrofcolumns);
%   transxarr2=zeros(nrofslices,nrofrows,nrofcolumns);
%   transyarr2=zeros(nrofslices,nrofrows,nrofcolumns);
%   problemflag=zeros(nrofslices,nrofrows,nrofcolumns);
%   
%   for row=1:1:nrofrows
%     for column=1:1:nrofcolumns
%       firstslice=1;
%       %for slice=1:1:nrofslices-1
%       for slice=startslicenr+1:1:startslicenr+nrofslices-1
%         slice
%         
%         if firstslice==1
%           firstslice=0;
%           name=sprintf(basename,slice-1); filename=sprintf('%s%s',indirectory,name);
%           img1=imread(filename);
%           %name=sprintf('../../SMG/scaled/SMG_scaled8_%d_R%dC%d.png',slice+startslicenr-1,row,column);
%           %img1=imread(name);
%           %img1=bandpass2(img1,0,100);
%           msimg1=downscalemean(img1,scale)/255;
%           msimg1=msimg1-mean(mean(msimg1));
%         else
%           img1=img2;
%           msimg1=msimg2;
%         end;
%         
%         name=sprintf(basename,slice); filename=sprintf('%s%s',indirectory,name);
%         img2=imread(filename);
%         %name=sprintf('../../SMG/scaled/SMG_scaled8_%d_R%dC%d.png',slice+startslicenr,row,column);
%         %img2=imread(name);
%         %img2=bandpass2(img2,0,100);
%         msimg2=downscalemean(img2,scale)/255;
%         msimg2=msimg2-mean(mean(msimg2));
%         
%         startangle=rotarr(slice,row,column)-5;
%         endangle=rotarr(slice,row,column)+5;
%         nrofangles=100;
%         rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
%         
%         figure(5);
%         imagesc(msimg1);
%         axis square;
%         colormap('gray');
%         figure(6);
%         imagesc(msimg2);
%         axis square;
%         colormap('gray');
%         
%         [transx,transy,corrval]=getrottransarrays(msimg1,msimg2,rots,'linear',0);
%         [C,D]=max(corrval);
%         rotarr2(slice,row,column)=rots(D);
%         if rotarr2(slice,row,column)>180
%           rotarr2(slice,row,column)=rotarr(slice,row,column)-360;
%         end;
%         %corrarr2(slice,:)=corrval;
%         transxarr2(slice,row,column)=(transx(D)+0.5)*scale;
%         transyarr2(slice,row,column)=(transy(D)+0.5)*scale;
%         
%         figure(9);
%         plot(rots,corrval);
%         grid on;
%         txt=sprintf('Slice %d compared to Slice %d',slice,slice-1);
%         title(txt);
%         xlabel('Rotation angle');
%         ylabel('Maximal correlation');
%         
%         mc=max(corrval);
%         if (corrval(1)==mc) || (corrval(end)==mc)
%           problemflag(slice)=1;
%         end;
%         
%         %       s=floor(size(msimg2,1)*sqrt(2)/2);
%         %       %visualize
%         %       figure(1);
%         %       imagesc(msimg1);
%         %       colormap(gray);
%         %       axis square;
%         %
%         %       figure(2);
%         %       rmsimg2=rotateimage(msimg2,rotarr2(slice,row,column)*pi/180,s,s,'linear');
%         %       imagesc(rmsimg2);
%         %       colormap(gray);
%         %       axis square;
%         
%         s=floor(size(img2,1)*sqrt(2)/2);
%         figure(3);
%         x1=(size(img1,2)/2)-(s/2)+transxarr2(slice,row,column); %transx is horizontal, transy is vertical.
%         y1=(size(img1,1)/2)-(s/2)+transyarr2(slice,row,column);
%         x2=min(size(img1,2),x1+s);
%         y2=min(size(img1,1),y1+s);
%         x1=max(x1,1); y1=max(y1,1); %lower image clipping
%         imagesc(img1(y1:y2,x1:x2));
%         axis([1 s 1 s]);
%         colormap(gray);
%         axis square;
%         
%         figure(4);
%         rimg2=rotateimage(double(img2),rotarr2(slice,row,column)*pi/180,s,s,'linear');
%         imagesc(rimg2);
%         colormap(gray);
%         axis square;
%         
%         pause(1);
%         %input('Press return...');
%       end;
%     end;
%   end;
%   
%   save('ROI_rigid3.mat','rotarr2','transxarr2','transyarr2','problemflag');
% else
%   load('ROI_rigid3.mat');
% end;
% 
% for row=1:1:nrofrows
%   for column=1:1:nrofcolumns
%     figure(10+(row-1)*nrofcolumns+(column-1));
%     subplot(3,1,1);
%     hold on;
%     plot(squeeze(rotarr2(:,row,column)),'g');
%     grid on;
%     title('Relative Rotation');
%     xlabel('Slice');
%     hold off;
%     
%     subplot(3,1,2);
%     hold on;
%     plot(squeeze(transxarr2(:,row,column)*scale),'g');
%     grid on;
%     title('Relative X Translation');
%     xlabel('Slice');
%     hold off;
%     
%     subplot(3,1,3);
%     hold on;
%     plot(squeeze(transyarr2(:,row,column)*scale),'g');
%     grid on;
%     title('Relative Y Translation');
%     xlabel('Slice');
%     hold off;
%   end;
% end;
% 
% %display progression of rotations, derived from absolute and relative measures
% load('../ROI_absrot_fft.mat');
% figure(19);
% plot(absrotang);
% hold on;
% plot(cumsum(-rotarr2)+absrotang(1),'r');
% grid on;
% hold off;
% 
% 
% 
% %compute relative rigid translation/rotation transformation matrices
% relrigid=zeros(nrofslices,3,3);
% for s=1:1:nrofslices
%   relrigid(s,1,1)=cos(rotarr2(s)*pi/180);
%   relrigid(s,1,2)=-sin(rotarr2(s)*pi/180);
%   relrigid(s,2,1)=sin(rotarr2(s)*pi/180);
%   relrigid(s,2,2)=cos(rotarr2(s)*pi/180);
%   relrigid(s,1,3)=transxarr2(s);
%   relrigid(s,2,3)=transyarr2(s);
%   relrigid(s,3,3)=1;
% end;
% 
% %compute absolute rigid translation/rotation transformation matrices
% absrigid=zeros(nrofslices,3,3);
% absrigid(1,:,:)=relrigid(1,:,:);
% for s=2:1:nrofslices
%   absrigid(s,:,:)=squeeze(absrigid(s-1,:,:))*squeeze(relrigid(s,:,:));
% end;
% 
% save('ROI_rigid3.mat','rotarr2','transxarr2','transyarr2','problemflag','relrigid','absrigid');
% 
% %compute overall bounding box assuming images of 1024^2
% for s=1:1:nrofslices
%   sbbox=computeboundingbox(1,1,1024,1024,squeeze(absrigid(s,:,:)));
%   if s==1
%     bbox=sbbox;
%   else
%     bbox(1)=min([bbox(1) sbbox(1)]);
%     bbox(2)=min([bbox(2) sbbox(2)]);
%     bbox(3)=max([bbox(3) sbbox(3)]);
%     bbox(4)=max([bbox(4) sbbox(4)]);
%   end;
% end;
% 
% 
% %%%%%%% Display result
% scale=16;
% for row=1:1:nrofrows
%   for column=1:1:nrofcolumns
%     firstslice=1;
%     %for slice=1:1:nrofslices-1
%     for slice=startslicenr+1:1:startslicenr+nrofslices-1
%       slice
%       
%       if firstslice==1
%         firstslice=0;
%         name=sprintf(basename,slice-1); filename=sprintf('%s%s',indirectory,name);
%         img1=imread(filename);
%         msimg1=downscalemean(img1,scale)/255;
%         msimg1=msimg1-mean(mean(msimg1));
%       else
%         img1=img2;
%         msimg1=msimg2;
%       end;
%       
%       name=sprintf(basename,slice); filename=sprintf('%s%s',indirectory,name);
%       img2=imread(filename);
%       msimg2=downscalemean(img2,scale)/255;
%       msimg2=msimg2-mean(mean(msimg2));
%       
%       s=floor(size(img2,1)*sqrt(2)/2);
%       figure(7);
%       x1=floor((size(img1,2)/2)-(s/2)+(transxarr2(slice,row,column))/scale); %transx is horizontal, transy is vertical.
%       y1=floor((size(img1,1)/2)-(s/2)+(transyarr2(slice,row,column))/scale);
%       x2=min(size(img1,2),x1+s);
%       y2=min(size(img1,1),y1+s);
%       xvofs=0; yvofs=0;
%       if (x1<1)
%         xvofs=x1;
%       end;
%       if (y1<1)
%         yvofs=y1;
%       end;
%       x1=max(x1,1); y1=max(y1,1); %upper-left image clipping
%       imagesc(img1(y1:y2,x1:x2));
%       axis([xvofs+1 xvofs+s yvofs+1 yvofs+s]);
%       colormap(gray);
%       axis square;
%       txt=sprintf('Slice %d Row %d Column %d',slice-1, row, column);
%       title(txt);
%       
%       figure(8);
%       rimg2=rotateimage(double(img2),rotarr2(slice,row,column)*pi/180,s,s,'linear');
%       imagesc(rimg2);
%       colormap(gray);
%       axis square;
%       txt=sprintf('Slice %d Row %d Column %d',slice, row, column);
%       title(txt);
%       
%       %pause(1);
%       input('Press return...');
%     end;
%   end;
% end;
