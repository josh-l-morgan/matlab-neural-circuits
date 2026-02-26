function param=pipeline_correctabsrotbetween(param)
%A function that detects outliers of estimated orientation and attempts to correct them
%Outliers are detected with respect to the average absolute orientation
%For use with the pipeline
%By Daniel Berger, June 15th 2009

lparam=param.correctabsrotbetween;

absrotang=param.absrot.uncorrang;
absrotpwr=param.absrot.uncorrpwr;
absrotang2=absrotang;
absrotpwr2=absrotpwr;
nooutliermean=zeros(param.nrofrows,param.nrofcolumns);
nooutlierstd=zeros(param.nrofrows,param.nrofcolumns);

disp('Absolute rotation Outlier detection...');
for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    txt=sprintf('Row %d, Column %d:',row,column);
    disp(txt);
    rcabsrotang=squeeze(absrotang(:,row,column));
    outlierspos=abs(rcabsrotang-mean(rcabsrotang))>std(rcabsrotang)*lparam.outlierthreshold;
    nrofoutliers=sum(outlierspos);
    if (nrofoutliers==0)
      disp('  No outliers found.');
      nooutliermean(row,column)=mean(rcabsrotang(outlierspos==0));
      nooutlierstd(row,column)=std(rcabsrotang(outlierspos==0));
    else
      disp('  Outliers found at:');
      disp(find(outlierspos));
      disp('  Refining outlier detection...');
      nooutliermean(row,column)=mean(rcabsrotang(outlierspos==0));
      nooutlierstd(row,column)=std(rcabsrotang(outlierspos==0));
      outlierspos2=abs(rcabsrotang-nooutliermean(row,column))>nooutlierstd(row,column)*lparam.outlierthreshold;
      nrofoutliers2=sum(outlierspos2);
      if (nrofoutliers2==0)
        disp('  No outliers found.');
      else
        disp('  Outliers found at:');
        disp(find(outlierspos2));
        
        disp('Computing local maxima around overall mean for outliers ...');
        seekmin=floor((nooutliermean(row,column)-nooutlierstd(row,column)*lparam.outlierthreshold)*10);
        seekmax=floor((nooutliermean(row,column)+nooutlierstd(row,column)*lparam.outlierthreshold)*10+1);
        
        for slice=1:1:param.nrofslices
          if outlierspos2(slice) %correct only outliers
            
            if (param.nrofrows==1)&&(param.nrofcolumns==1)
              name=sprintf(param.basescaledname,slice);
            else
              name=sprintf(param.basescaledname,slice,row,column);
            end;
            filename=sprintf('%s%s',param.scaleddir,name);
            %name=sprintf(param.basescaledname,slice); filename=sprintf('%s%s',param.scaleddir,name);
            %name=sprintf('../scaled/w7_ROIb_scaled16_s%02d.png',slice);
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              im=imread(filename);
              disp('Computing ...');
              
              if (size(im,1)~=size(im,2))||(mod(size(im,1),2)>0.5) %make divisible-by-two square image
                %(Modulo value should be 0 or 1, thresholding it at 0.5 makes this rounding-error resistent)
                dsize=min(size(im));
                if (mod(dsize,2)>0.5)
                  dsize=dsize-1;
                end;
              else
                dsize=size(im,1);
              end;
              cim=im(1:dsize,1:dsize); %crop image.
              fim=abs(fftshift(fft2(cim))); %do fft
              figure(2);
              imagesc(log(fim));
              hdsize=floor(dsize/2); %half of dsize
              lfreq=floor(lparam.lfreq*hdsize); %lowest frequency to consider
              hfreq=floor(lparam.hfreq*hdsize); %highest frequency to consider
              %make polar->cartesian transform using 1/10th degree resolution
              pim=topolar4(fim(2:dsize,2:dsize),([lfreq:hfreq])',([0:3599]*0.1*pi/180)','linear');
              figure(3);
              imagesc(log(pim));
              
              pim2=sum(pim); %integral over all relevant frequencies for each orientation
              pim2=pim2(1:1800)+pim2(1801:3600); %fold 0..360deg to 0..180deg)
              figure(5);
              plot([0:1799]/10,pim2);
              txt=sprintf('Slice %d',slice);
              title(txt);
              hold on;
              plot([seekmin/10 seekmin/10],[min(pim2) max(pim2)],'r');
              plot([seekmax/10 seekmax/10],[min(pim2) max(pim2)],'r');
              grid on;
              hold off;
              
              [pwr,ang]=max(pim2(seekmin:seekmax));
              absrotpwr2(slice,row,column)=pwr;
              absrotang2(slice,row,column)=(seekmin+ang)/10;
              %if lparam.stopeach
              input('Press RETURN to continue');
              %end;
              
%               tim=toc;
%               count=count+1;
%               average_time=tim/count
            else
              txt=sprintf('Warning: File %s not found.', filename);
              disp(txt);
            end;
          end;
        end;
      end;
    end;
  end;
end;
  
  param.correctabsrotbetween.ang=absrotang2;
  param.correctabsrotbetween.pwr=absrotpwr2;
  param.correctabsrotbetween.nooutliermean=nooutliermean;
  param.correctabsrotbetween.nooutlierstd=nooutlierstd;


figure(19);
count=1;
for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    subplot(param.nrofrows,param.nrofcolumns,count);
    plot(squeeze(absrotang(:,row,column)));
    hold on;
    plot(squeeze(absrotang2(:,row,column)),'g');
    grid on;
    t=sprintf('Uncorrected (blue) and corrected (green) orientation estimates (fft-based) R%dC%d',row,column);
    title(t);
    hold off;
    count=count+1;
  end;
end;