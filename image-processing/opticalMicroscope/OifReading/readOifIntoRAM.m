clear all

CPN = 'H:\LxA_GadGFP_Albino\Mouse2\B3\20X_2d.oif.files\'

stkChan = {[3 2 1] [5 4  2]};
useSlices = [];%[1:34];


mkdir(TPN)

dCPN = dir([CPN '*.tif'])
inams = {dCPN.name};
clear color slice
for i  = 1:length(inams)
    
    nam = inams{i};
    
    c = regexp(nam,'C');
    d = regexp(nam,'.tif');
    z = regexp(nam,'Z');
    
    if isempty(z)
            color(i) = str2num(nam(c+1:d-1));
            slice(i) = 1;
    else
    color(i) = str2num(nam(c+1:z-1));
    slice(i) = str2num(nam(z+1:d-1));
    end
    
    
end




if isempty(useSlices)
    slices = unique(slice);
else
    slices = useSlices;
end



if 0
    for s = 1:length(stkChan)
        TPN_I = sprintf('%sI_%d\\',TPN,s);
        if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
        mkdir(TPN_I)
        for i = 1:length(slices)
            
            clear Ic
            targ = find((slice == slices(i)) & ( color == stkChan{s}(1)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = medfilt2(I,[3 3]);
                I = I/256;
                Ic(:,:,1) = uint8(I);
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(2)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = I/256;
                I = medfilt2(I,[3 3]);
                Ic(:,:,2) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(3)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = medfilt2(I,[3 3]);
                I = I/256;
                Ic(:,:,3) = I;
            end
            if size(Ic,3)<3
                Ic(1,1,3) = 0;
            end
            
            image(uint8(Ic*100))
            i
            pause
            pause(.01)
            imwrite(uint8(Ic),sprintf('%s%05.0f.tif',TPN_I,slices(i)))
            
        end
    end
    
end


if 1
    
    testI = imread([CPN inams{1}]);
    [ys xs zs] = size(testI);
    for s = 1:length(stkChan)
        TPN_I = sprintf('%sI_%d\\',TPN,s);
        if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
        mkdir(TPN_I)
        Ic1 = zeros(ys,xs,length(slices),'double');
        Ic2 = Ic1; Ic3 = Ic1;
        for i = 1:length(slices)
            
            targ = find((slice == slices(i)) & ( color == stkChan{s}(1)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                I = medfilt2(I,[3 3]);
                Ic1(:,:,i) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(2)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                I = medfilt2(I,[3 3]);
                Ic2(:,:,i) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(3)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                I = medfilt2(I,[3 3]);
                Ic3(:,:,i) = I;
            end
        end
        
        Ic1 = Ic1 * 256/max(Ic1(:));
        Ic2 = Ic2 * 256/max(Ic2(:));
        Ic3 = Ic3 * 256/max(Ic3(:));
        
        Imax1 = max(Ic1,[],3);
        Imax2 = max(Ic2,[],3);
        Imax3 = max(Ic3,[],3);
        
        Imax = (cat(3,Imax1,Imax2,Imax3));
        
        Imean1 = mean(Ic1,3);
        Imean2 = mean(Ic2,3);
        Imean3 = mean(Ic3,3);
        
        Imean = cat(3,Imean1,Imean2,Imean3);
        
        imwrite(uint8(Imax),sprintf('%smax_stk%d.tif',TPN,s))
        
        
        
        
        
            TPN_I = sprintf('%sI_%d\\',TPN,s);
            if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
            mkdir(TPN_I)
            for i = 1:length(slices)
                Ic = cat(3,Ic1(:,:,i),Ic2(:,:,i),Ic3(:,:,i));
                image(uint8(Ic*100))
                pause(.01)
                imwrite(uint8(Ic),sprintf('%s%05.0f.tif',TPN_I,slices(i)))
            end
        
        
    end
    
    
end












