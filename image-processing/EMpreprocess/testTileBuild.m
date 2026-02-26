
%% Tile test reader

clear all
SPNsec = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2_reDice7\303\'

for s = 0:8;
    
    SPN = sprintf('%s%d\\',SPNsec,s);
    
    dSPN = dir([SPN '*.png']);
    inams = {dSPN.name};
    
    
    
    dI = dir([SPN '*.png']);
    inams = {dI.name};
    rs = zeros(length(inams),1);
    cs = rs;
    Ir = {};
    
    for i = 1:length(inams)
        
        nam = inams{i};
        u = regexp(nam,'_');
        d = regexp(nam,'.png');
        rs(i) = str2num(nam(1:u-1));
        cs(i) = str2num(nam(u+1:d-1));
        Ir{i} = imread([SPN nam]);
    end
    
    
    
    [ys xs] = size(Ir{1});
    maxR = max(rs)
    maxC = max(cs)
    
    I = zeros(maxR*ys,maxC * xs);
    size(I)
    
    colormap gray(256)
    
    if 0
        for i = 1:length(Ir)
            Is = Ir{i};
            if ~isempty(Is)
                image(Is)
                pause
            end
        end
    end
    
    
    
    for i = 1:length(Ir)
        disp(sprintf('showing %d of %d',i,length(Ir)))
        
        [ys xs] = size(Ir{i});
        
        
        startY = (rs(i) - 1) * ys + 1;
        stopY =  rs(i) * ys;
        
        startX = (cs(i) - 1) * xs + 1;
        stopX = cs(i) * xs;
        
        I(startY:stopY, startX:stopX) = Ir{i};
        
    end
    
    Ia{s+1} = I;
    
    imshow(uint8(I))
    pause(.1)
end


% %%%%%%%
% Ic = imresize(Ia{1},1/4);
% Ic(:,:,2) = imresize(Ia{2},1/2);
% Ic(:,:,3) = imresize(Ia{3},1);
% image(Ic)

subplot(2,2,1)
image(Ia{2});
subplot(2,2,2)
image(Ia{3});

subplot(2,2,3)
image(Ia{4});
subplot(2,2,4)
image(Ia{5});






