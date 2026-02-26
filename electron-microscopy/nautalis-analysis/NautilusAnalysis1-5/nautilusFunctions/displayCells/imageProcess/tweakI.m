function[I] = tweakI(I,viewProps);  
%% general image contrast tweaker


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



    if isfield(viewProps,'keepRat')
        keepRat = viewProps.keepRat;
        I = uint8keepRat(I)*keepRat + uint8(I)*(1-keepRat);
    end





