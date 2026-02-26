function[mfIs, listYX] = funSubFFT(Ir);
    
    sS = 100;

    if size(Ir,1)>500 & size(Ir,2)>500;
        I = Ir(100:end-100,10:end-100);
    else
        I = Ir;
    end
    [ys xs] = size(I);
    %%Fix Contrast
    I = I - min(I(:));
    I = I * 255/double(max(I(:)));
    
    yStep = 1:sS:ys-sS;
    xStep = 1:sS:xs-sS;
    fftMos = zeros(length(yStep),length(xStep));
    mfIs = zeros(length(yStep) * length(xStep),round(sS/2) - 1);
    c = 0;
    for y = 1:length(yStep)
        for x = 1:length(xStep)
            c = c+1;
            subI = I(yStep(y):yStep(y)+sS-1,xStep(x):xStep(x)+sS-1);
%             subI = subI - min(subI(:));
%             subI = subI * 255/max(subI(:));
            subI = subI - mean(subI(:));
            fftI = abs(fft(subI,[],1));
            
            mfIa = mean(fftI,2);
            mfI = mfIa(2:fix(length(mfIa)/2));
            bgFFT = mean(mfI(end-10:end));
            mfI = mfI - bgFFT;
            mfIs(c,:) = mfI;
            
            
            listYX(c,:) = [y x];
            
        end
    end
    
    
    

