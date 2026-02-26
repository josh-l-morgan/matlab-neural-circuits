  function param=pipeline_checkrelrotconsistency(param)
%A function to check whether the relative orientations of consecutive slices

lparam=param.checkrelrotconsistency;

if (param.nrofrows>1) || (param.nrofcolumns>1)
  
  rotarr=param.relrot.rotarr1;
  transxarr=param.relrot.transxarr1;
  transyarr=param.relrot.transyarr1;
  
  startangle=0; endangle=360; nrofangles=lparam.nrofangles;
  rots=startangle+[0:1:nrofangles-1]*(endangle-startangle)/nrofangles;
  
  stdrotarr1=zeros(param.nrofslices,1);
  for slice=1:1:param.nrofslices
    val=reshape(param.relrot.rotarr1(slice,:,:),1,param.nrofrows*param.nrofcolumns);
    stdrotarr1(slice)=std(val);
  end;
  
  figure(14);
  p1=plot(stdrotarr1);
  grid on;
  title('Relative rotations consistency check');
  xlabel('Slice No.');
  ylabel('STD of relative rotations for set of tiles in each slice');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Re-estimate the relative orientation for all slices in which a high STD is detected
  for slice=1:1:param.nrofslices
    if stdrotarr1(slice)>lparam.reestimatestdthres
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns
          %Load previous-slice tile
          name=sprintf(param.basescaledname,slice-1,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img1=imread(filename);
          if lparam.dofiltering
            img1=bandpass2(img1,lparam.lowestpixpercyc,lparam.highestpixpercyc);
          end;
          msimg1=downscalemean(img1,lparam.scale)/255;
          msimg1=msimg1-mean(mean(msimg1));
          
          %Load this-slice tile
          name=sprintf(param.basescaledname,slice,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);          
          disp(filename);
          img2=imread(filename);
          if lparam.dofiltering
            img2=bandpass2(img2,lparam.lowestpixpercyc,lparam.highestpixpercyc);
          end;
          msimg2=downscalemean(img2,lparam.scale)/255;
          msimg2=msimg2-mean(mean(msimg2));

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
          figure(2);
          imagesc(msimg2);
          colormap(gray);
          axis square;
          
          s=floor(size(img2,1)*sqrt(2)/2);
          figure(3);
          x1=floor((size(img1,2)/2)-(s/2)+transxarr(slice)); %transx is horizontal, transy is vertical.
          y1=floor((size(img1,1)/2)-(s/2)+transyarr(slice));
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
      end;
    end;
  end;
  
  param.checkrelrotconsistency.rotarr1=rotarr;
  param.checkrelrotconsistency.transxarr1=transxarr;
  param.checkrelrotconsistency.transyarr1=transyarr;
  
  stdrotarr2=zeros(param.nrofslices,1);
  for slice=1:1:param.nrofslices
    val=reshape(rotarr(slice,:,:),1,param.nrofrows*param.nrofcolumns);
    stdrotarr2(slice)=std(val);
  end;
  
  figure(14);
  hold on;
  p2=plot(stdrotarr2,'g');
  hold off;
  
  figure(15);
  plot(reshape(rotarr,param.nrofslices,param.nrofrows*param.nrofcolumns));
  grid on;
  title('Individual rotation estimates of the slice tiles');
  xlabel('Slice No.');
  ylabel('Rotation relative to previous slice [deg]');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Check if there are still outliers; If yes, use Median; 
  %(if tie, use smaller one <-- not implemented at the moment)
  
  for slice=1:1:param.nrofslices
    if stdrotarr2(slice)>lparam.reestimatestdthres
      rot=median(reshape(rotarr(slice,:,:),1,size(rotarr,2)*size(rotarr,3)));
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns
          rotarr(slice,row,column)=rot;
          
          %re-estimate translations, using the given rotation
          
          %Load previous-slice tile
          name=sprintf(param.basescaledname,slice-1,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);
          disp(filename);
          img1=imread(filename);
          if lparam.dofiltering
            img1=bandpass2(img1,lparam.lowestpixpercyc,lparam.highestpixpercyc);
          end;
          msimg1=downscalemean(img1,lparam.scale)/255;
          msimg1=msimg1-mean(mean(msimg1));
          
          %Load this-slice tile
          name=sprintf(param.basescaledname,slice,row,column);
          filename=sprintf('%s%s',param.scaleddir,name);          
          disp(filename);
          img2=imread(filename);
          if lparam.dofiltering
            img2=bandpass2(img2,lparam.lowestpixpercyc,lparam.highestpixpercyc);
          end;
          msimg2=downscalemean(img2,lparam.scale)/255;
          msimg2=msimg2-mean(mean(msimg2));

          txt=sprintf('Computing maximal cross-correlation peak for relative rotation of %d degrees...',rot);
          disp(txt);
          [transx,transy,corrval]=getrottransarrays(msimg1,msimg2,rot,'linear',0);
          [C,D]=max(corrval);
          rotarr(slice,row,column)=rot; %rots(D);
          if rotarr(slice,row,column)>180
            rotarr(slice,row,column)=rotarr(slice,row,column)-360;
          end;
          %corrarr(slice,:)=corrval;
          transxarr(slice,row,column)=(transx+0.5)*lparam.scale; %(transx(D)+0.5)*lparam.scale;
          transyarr(slice,row,column)=(transy+0.5)*lparam.scale; %(transy(D)+0.5)*lparam.scale;
        end;
      end;
    end;
  end;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % rotarr, transxarr and transyarr are the final corrected coarse estimates.
  param.checkrelrotconsistency.rotarr=rotarr;
  param.checkrelrotconsistency.transxarr=transxarr;
  param.checkrelrotconsistency.transyarr=transyarr;
 
  stdrotarr3=zeros(param.nrofslices,1);
  for slice=1:1:param.nrofslices
    val=reshape(rotarr(slice,:,:),1,param.nrofrows*param.nrofcolumns);
    stdrotarr3(slice)=std(val);
  end;
  
  figure(14);
  hold on;
  p3=plot(stdrotarr3,'r');
  hold off;
  legend([p1 p2 p3],'Original estimates','After outlier correction','After median correction');
  
else
  disp('Error: There is only one row and one column, so no comparison is possible.');
end;