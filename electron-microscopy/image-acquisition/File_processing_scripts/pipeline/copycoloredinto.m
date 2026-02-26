function targetimage=copycoloredinto(targetimage,source,channels)
%A function which copies the one-channel image source in the given channels of the targetimage
%and returns the targetimage. Channels is three numbers, each 1or 0.
%This is not a memory-efficient version since it copies the target image.
%targetimage and source need to have the same (x,y) size.
%source is 2d, targetimage 3d (2d + 3 channels)
%By Daniel Berger for MIT-BCS Seung, June 12th 2009

if (size(targetimage,1)==size(source,1)) && (size(targetimage,2)==size(source,2))
  for channel=1:1:3
    if (channels(channel))>0
      for row=1:1:size(targetimage,1)
        for column=1:1:size(targetimage,2)
          if isfinite(source(row,column))
            targetimage(row,column,channel)=source(row,column);
          end;
        end;
      end;
    end;
  end;
else
  disp('copycoloredinto.m ERROR: source and target have different sizes.');
end;