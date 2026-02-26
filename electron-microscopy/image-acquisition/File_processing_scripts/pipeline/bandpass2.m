function fimage = bandpass2(image, lowcut, highcut)
%Does a bandpass filtering using a ring mask on a 2d-fft of the image.
%lowcut and highcut are the cutoff frequencies of the bandpass filter,
%in pixels/cycle (lowcut is the lowest pix/cycle, so the highest freq)
%By Daniel Berger, March 11 2009 (MIT-BCS Seung)

lowcutsq=lowcut*lowcut;
highcutsq=highcut*highcut;

fftimg=fftshift(fft2(image));

%lowest freq is at size(image)/2+1
ringmask=zeros(size(image));

%coordinates of constant value component in the fft image
midy=floor(size(image,1)/2+1);
midx=floor(size(image,2)/2+1);

%compute ring mask
for y=1:1:size(ringmask,1)
  dy=abs(midy-y);
  if dy==0
    for x=1:1:size(ringmask,2)
      dx=abs(midx-x);
      pixpp_x=size(ringmask,2)/dx;
      pixpp_xsq=pixpp_x*pixpp_x;
      if (pixpp_xsq>=lowcutsq) && (pixpp_xsq<=highcutsq)
        ringmask(y,x)=1;
      end;
    end;
  else
    pixpp_y=size(ringmask,1)/dy;
    pixpp_ysq=pixpp_y*pixpp_y;
    if (pixpp_ysq>=lowcutsq)
      for x=1:1:size(ringmask,2)
        dx=abs(midx-x);
        if dx==0
          if (pixpp_ysq>=lowcutsq) && (pixpp_ysq<=highcutsq)
            ringmask(y,x)=1;
          end;
        else
          pixpp_x=size(ringmask,2)/dx;
          pixpp_xsq=pixpp_x*pixpp_x;
          pixpp_xysq=(pixpp_xsq*pixpp_ysq)/(pixpp_x*pixpp_x+pixpp_y*pixpp_y);
          if (pixpp_xysq>=lowcutsq) && (pixpp_xysq<=highcutsq)
            ringmask(y,x)=1;
          end;
        end;
      end;
    end;
  end;
end;

%apply ring mask and transform back
fftimg=fftimg.*ringmask;
fimage = ifft2(ifftshift(fftimg));