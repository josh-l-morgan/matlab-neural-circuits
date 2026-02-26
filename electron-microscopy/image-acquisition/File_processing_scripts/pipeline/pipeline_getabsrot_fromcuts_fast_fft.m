%function [absrotang,absrotpwr]=pipeline_getabsrot_fromcuts_fast_fft(param)
function param=pipeline_getabsrot_fromcuts_fast_fft(param)
%A function to compute fast absolute image orientation estimates
%from knife cutting artifacts.
%This version is for usage in the pipeline.
%By Daniel Berger for MIT-BCS Seung, April 26 2009

lparam=param.absrot;


absrotang=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);
absrotpwr=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns);

tic;
count=0;

for slice=1:1:param.nrofslices
  for row=1:1:param.nrofrows
    for column=1:1:param.nrofcolumns
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
        grid on;
      
        [pwr,ang]=max(pim2);
        absrotpwr(slice,row,column)=pwr;
        absrotang(slice,row,column)=ang/10;
        if lparam.stopeach
          input('Press RETURN to continue');
        end;
        
        tim=toc;
        count=count+1;
        average_time=tim/count
      else
        txt=sprintf('Warning: File %s not found.', filename);
        disp(txt);
      end;
    end;
  end;
end;


figure(19);
count=1;
for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    subplot(param.nrofrows,param.nrofcolumns,count);
    plot(squeeze(absrotang(:,row,column)));
    hold on;
    grid on;
    t=sprintf('Uncorrected orientation estimates (fft-based) R%dC%d',row,column);
    title(t);
    xlabel('Slice No.');
    ylabel('Absolute orientation [deg]');
    %axis([0 60 0 180]);
    hold off;
    count=count+1;
  end;
end;

param.absrot.uncorrang=absrotang;
param.absrot.uncorrpwr=absrotpwr;
