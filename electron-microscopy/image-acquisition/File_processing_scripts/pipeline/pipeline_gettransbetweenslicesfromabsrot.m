function param=pipeline_gettransbetweenslicesfromabsrot(param)
%A function to compute relative translations of tiles between slices
%by using the results of the absolute orientation estimation
%To be used within the pipeline
%By Daniel Berger for MIT-BCS Seung, June 4th 2009

lparam=param.abstransbetweenslices;
absrotang=lparam.in_absrotang;

param.abstransbetweenslices.transx=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
param.abstransbetweenslices.transy=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
param.abstransbetweenslices.corrval=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);

for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    firstslice=1;
    for slice=2:1:param.nrofslices
      %//Load images.
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
        img1=single(img1)/255;
        if lparam.dofiltering
          img1=bandpass2(img1,lparam.lowestpixpercyc,lparam.highestpixpercyc);
        end;
      else
        img1=img2;
      end;
      
      if (param.nrofrows==1)&&(param.nrofcolumns==1)
        name=sprintf(param.basescaledname,slice);
      else
        name=sprintf(param.basescaledname,slice,row,column);
      end;
      filename=sprintf('%s%s',param.scaleddir,name);
      disp(filename);
      img2=imread(filename);
      img2=single(img2)/255;
      if lparam.dofiltering
        img2=bandpass2(img2,lparam.lowestpixpercyc,lparam.highestpixpercyc);
      end;
      
      %Rotate second image so that it is rotationally aligned with the first (to mimic relative orientation estimates)
      %[transx,transy,corrval]=getrottransarrays(img1,img2,(param.absrot.ang(slice)-param.absrot.ang(slice-1)),'linear');
      relang=absrotang(slice,row,column)-absrotang(slice-1,row,column);
      if (slice==305)
        [transx,transy,corrval]=getrottransarrays(img1,img2,relang,'linear',1);
        input('Press return...');
      else
        [transx,transy,corrval]=getrottransarrays(img1,img2,relang,'linear',0);
      end;
%      [C,D]=max(corrval);
%       rotarr(slice,row,column)=rots(D);
%       if rotarr(slice,row,column)>180
%         rotarr(slice,row,column)=rotarr(slice,row,column)-360;
%       end;
      %corrarr(slice,:)=corrval;
      param.abstransbetweenslices.transx(slice,row,column)=transx; %(transx(D)+0.5)*lparam.scale;
      param.abstransbetweenslices.transy(slice,row,column)=transy; %(transy(D)+0.5)*lparam.scale;
      param.abstransbetweenslices.corrval(slice,row,column)=corrval;
      
%       if slice==305
%         %visualize
% %         figure(1);
% %         imagesc(msimg1);
% %         colormap(gray);
% %         axis square;
% %         
% %         s=floor(size(msimg2,1)*sqrt(2)/2);
% %         rmsimg2=rotateimage(msimg2,rots(D)*pi/180,s,s,'linear');
% %         figure(2);
% %         imagesc(rmsimg2);
% %         colormap(gray);
% %         axis square;
%         
%         s=floor(size(img2,1)*sqrt(2)/2);
%         figure(3);
%         x1=(size(img1,2)/2)-(s/2)+transxarr(slice); %transx is horizontal, transy is vertical.
%         y1=(size(img1,1)/2)-(s/2)+transyarr(slice);
%         x2=min(size(img1,2),x1+s);
%         y2=min(size(img1,1),y1+s);
%         x1=max(x1,1); y1=max(y1,1); %lower image clipping
%         imagesc(img1(y1:y2,x1:x2));
%         axis([1 s 1 s]);
%         colormap(gray);
%         axis square;
%         txt=sprintf('Slice %d, which is compared with Slice %d',slice-1,slice);
%         title(txt);
%       
%         figure(4);
%         rimg2=rotateimage(double(img2),rots(D)*pi/180,s,s,'linear');
%         imagesc(rimg2);
%         colormap(gray);
%         axis square;
%         txt=sprintf('Slice %d, which is compared with Slice %d',slice,slice-1);
%         title(txt);
%         
%         %pause(1);
%         %if lparam.stopeach
%           input('Return..');
%         %end;
%       end;
      %Compute best-fitting translation (as in relative orientation/translation estimation)
    end;
  end;
end;


