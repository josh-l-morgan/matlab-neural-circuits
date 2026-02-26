function targetimage=addimageinto(targetimage,source)
%A function which adds the one-channel image source to the targetimage
%and returns the targetimage.
%This is not a memory-efficient version since it copies the target image.
%targetimage and source need to have the same (x,y) size.
%source and target are 2d, (1 channel)
%By Daniel Berger for MIT-BCS Seung, June 12th 2009

if (size(targetimage,1)==size(source,1)) && (size(targetimage,2)==size(source,2))
  for row=1:1:size(targetimage,1)
    for column=1:1:size(targetimage,2)
      if isfinite(source(row,column))
        targetimage(row,column)=targetimage(row,column)+source(row,column);
      end;
    end;
  end;
else
  disp('addimageinto.m ERROR: source and target have different sizes.');
end;