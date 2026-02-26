function param=pipeline_correctabsrotwithin(param)
%This function compares the absolute rotation estimates between tiles in the 
%same slice, assuming that the variability within a slice should be minimal.
%For use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 15 2009

lparam=param.correctabsrotwithin;
absrotang=lparam.in_absrotang;
absrotpwr=lparam.in_absrotpwr;

disp('Outlier detection...');
outlierspos=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
  
%Look for high variance of orientations within slice
for slice=1:1:param.nrofslices
  rcabsrotang=squeeze(absrotang(slice,:,:));
  rcara_mean=mean(mean(rcabsrotang));
  rcara_std=std(reshape(rcabsrotang,1,size(rcabsrotang,1)*size(rcabsrotang,2)));
  if lparam.recomputealltiles==1
    %If the STD is too high, mark all tiles for re-evaluation
    if rcara_std>lparam.outlierthreshold
      outlierspos(slice,:,:)=1;
    end;
  else
    %Identify which tiles are the outlier tiles
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        if abs(rcabsrotang(row,column)-rcara_mean)>lparam.outlierthreshold
          outlierspos(slice,row,column)=1;
        end;
      end;
    end;
  end;
end;

nrofoutliers=sum(sum(sum(outlierspos)));

figure(20);
plot(reshape(absrotang,size(absrotang,1),size(absrotang,2)*size(absrotang,3)));
hold on;
for slice=1:1:param.nrofslices
  for row=1:1:param.nrofrows
    for column=1:1:param.nrofcolumns
      if outlierspos(slice,row,column)>0
        plot(slice,absrotang(slice,row,column),'ko');
      end;
    end;
  end;
end;
hold off;
grid on;
xlabel('Slice No.');
ylabel('Absolute rotation estimate (deg.)');
txt=sprintf('Showing %d outliers',nrofoutliers);
title(txt);

if (nrofoutliers==0)
  disp('  No outliers found.');
  nooutliermean=mean(absrotang(outlierspos==0));
  nooutlierstd=std(absrotang(outlierspos==0));
else
  disp('  Outliers found at:');
  disp(find(outlierspos));
  
  %Re-compute absolute orientation of tiles with outliers
  %Search range is tile mean+/- 
  
  for slice=1:1:param.nrofslices
    rcabsrotang=squeeze(absrotang(slice,:,:));
    rcara_mean=mean(mean(rcabsrotang));
    rcara_std=std(reshape(rcabsrotang,1,size(rcabsrotang,1)*size(rcabsrotang,2)));
    startang=floor((rcara_mean-rcara_std*lparam.outlierthreshold));
    endang=floor((rcara_mean+rcara_std*lparam.outlierthreshold)+1);
    
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        if outlierspos(slice,row,column)>0
          %Re-estimate this tile.
          if (param.nrofrows==1)&&(param.nrofcolumns==1)
            name=sprintf(param.baserawname,slice);
          else
            name=sprintf(param.baserawname,slice,row,column);
          end;
          filename=sprintf('%s%s',param.rawdir,name);
          
%           if (param.nrofrows==1)&&(param.nrofcolumns==1)
%             name=sprintf(param.basescaledname,slice);
%           else
%             name=sprintf(param.basescaledname,slice,row,column);
%           end;
%           filename=sprintf('%s%s',param.scaleddir,name);
          if exist(filename,'file')
            txt=sprintf('Loading image %s ...',filename);
            disp(txt);
            im=imread(filename);
            im=double(im)/255;
            im=imresize(im,1/4);
      
            [pwr,ang]=getrotation(im,0,0.1,180,15,21);
            %input('Press return...');
            
            disp('Computing local maxima around overall mean for outliers ...');
            %[pwr,ang]=getrotation(im,startang,0.01,endang,15,31);
            [pwr,ang]=getrotation(im,180-endang,0.01,180-startang,15,31);
            absrotpwr(slice,row,column)=pwr;
            absrotang(slice,row,column)=180-ang;
            if (lparam.stopeach==1)
              input('Press return...');
            end;
          end;
        end;
      end;
    end;
  end;
end;

param.correctabsrotwithin.absrotang=absrotang;
param.correctabsrotwithin.absrotpwr=absrotpwr;

%           disp('Computing local maxima around overall mean for outliers ...');
%           seekmin=floor((nooutliermean(row,column)-nooutlierstd(row,column)*lparam.outlierthreshold)*10);
%           seekmax=floor((nooutliermean(row,column)+nooutlierstd(row,column)*lparam.outlierthreshold)*10+1);
%           
%           for slice=1:1:param.nrofslices
%             if outlierspos2(slice) %correct only outliers
%               
%               if (param.nrofrows==1)&&(param.nrofcolumns==1)
%                 name=sprintf(param.basescaledname,slice);
%               else
%                 name=sprintf(param.basescaledname,slice,row,column);
%               end;
%               filename=sprintf('%s%s',param.scaleddir,name);
%               %name=sprintf(param.basescaledname,slice); filename=sprintf('%s%s',param.scaleddir,name);
%               %name=sprintf('../scaled/w7_ROIb_scaled16_s%02d.png',slice);
%               if exist(filename,'file')
%                 txt=sprintf('Loading image %s ...',filename);
%                 disp(txt);
%                 im=imread(filename);
%                 disp('Computing ...');
%                 
%                 if (size(im,1)~=size(im,2))||(mod(size(im,1),2)>0.5) %make divisible-by-two square image
%                   %(Modulo value should be 0 or 1, thresholding it at 0.5 makes this rounding-error resistent)
%                   dsize=min(size(im));
%                   if (mod(dsize,2)>0.5)
%                     dsize=dsize-1;
%                   end;
%                 else
%                   dsize=size(im,1);
%                 end;
%                 cim=im(1:dsize,1:dsize); %crop image.
%                 fim=abs(fftshift(fft2(cim))); %do fft
%                 figure(2);
%                 imagesc(log(fim));
%                 hdsize=floor(dsize/2); %half of dsize
%                 lfreq=floor(lparam.lfreq*hdsize); %lowest frequency to consider
%                 hfreq=floor(lparam.hfreq*hdsize); %highest frequency to consider
%                 %make polar->cartesian transform using 1/10th degree resolution
%                 pim=topolar4(fim(2:dsize,2:dsize),([lfreq:hfreq])',([0:3599]*0.1*pi/180)','linear');
%                 figure(3);
%                 imagesc(log(pim));
%                 
%                 pim2=sum(pim); %integral over all relevant frequencies for each orientation
%                 pim2=pim2(1:1800)+pim2(1801:3600); %fold 0..360deg to 0..180deg)
%                 figure(5);
%                 plot([0:1799]/10,pim2);
%                 txt=sprintf('Slice %d',slice);
%                 title(txt);
%                 hold on;
%                 plot([seekmin/10 seekmin/10],[min(pim2) max(pim2)],'r');
%                 plot([seekmax/10 seekmax/10],[min(pim2) max(pim2)],'r');
%                 grid on;
%                 hold off;
%                 
%                 [pwr,ang]=max(pim2(seekmin:seekmax));
%                 absrotpwr2(slice,row,column)=pwr;
%                 absrotang2(slice,row,column)=(seekmin+ang)/10;
%                 %if lparam.stopeach
%                 input('Press RETURN to continue');
%                 %end;
%                 
%                 tim=toc;
%                 count=count+1;
%                 average_time=tim/count
%               else
%                 txt=sprintf('Warning: File %s not found.', filename);
%                 disp(txt);
%               end;
%             end;
%           end;
%         end;
%       end;
%     end;
% end;
% 
% 
% 
