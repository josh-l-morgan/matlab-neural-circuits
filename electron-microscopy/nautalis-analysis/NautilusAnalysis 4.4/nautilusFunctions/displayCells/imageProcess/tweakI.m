function[I] = tweakI(Ir,viewProps);  
%% general image contrast tweaker

I = Ir;
I = double(I);
    if isfield(viewProps,'dilate')
        SE = strel('disk',viewProps.dilate);
        I = imdilate(I,SE);
    end


    if isfield(viewProps,'contrast')
        I = I * viewProps.contrast;
    end


    if isfield(viewProps,'gamma')
        gamma = viewProps.gamma;
        I = I.^gamma * 256/256^gamma;
    end



%     if isfield(viewProps,'keepRat')
%         keepRat = viewProps.keepRat;
%         clipped = I;
%         clipped(clipped>255) = 255;
%         I = (I*255/max(I(:)) * keepRat) + clipped * (1-keepRat);
%     end
% 
% image(I)
% image(Ir);
% 


