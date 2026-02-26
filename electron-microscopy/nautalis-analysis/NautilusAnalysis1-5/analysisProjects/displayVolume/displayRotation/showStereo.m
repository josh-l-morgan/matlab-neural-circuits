function[I_stereo] = showStereo(viewProps)

viewProps.degRot = -2;
I_stereo1 = stereoCellsAndMoreFull(viewProps);
% image(uint8(I_stereo1*3))
[y1 x1] = find(sum(I_stereo1,3));


viewProps.degRot = +2;
I_stereo2 = stereoCellsAndMoreFull(viewProps);
[y2 x2] = find(sum(I_stereo2,3));

clipY = [min([y1;y2]) max([y1;y2])];
I_stereo1 = I_stereo1(clipY(1):clipY(2),min(x1):max(x1),:);
I_stereo2 = I_stereo2(clipY(1):clipY(2),min(x2):max(x2),:);



I_stereo = cat(2,I_stereo1,I_stereo2);
% image(uint8(I_stereo))

I_stereo = cat(2,I_stereo2,I_stereo1);

% 
% Ist = uint8keepRat(I_stereo*.3) + uint8(I_stereo*.1);
% Ist = uint8keepRat(I_stereo*.03) + uint8(I_stereo*.01);
% 
% image(Ist*10);



