


SPN = 'E:\affFull_fullRes_mipd\';

secDir = dir(SPN);
secFold = {secDir([secDir.isdir]>0).name};

iSize = 1024;

for s = 1:length(secFold)
    sec = str2num(secFold{s})
    if ~isempty(sec)
        mipDir = [SPN secFold{s} '\'];
        inDir = dir(mipDir);
        mipLevs = {inDir.name};
        maxLev = 4;
%         for i = 1:length(mipLevs);
%             if ~isempty(str2num(mipLevs{i}));
%                 maxLev = max(maxLev,str2num(mipLevs{i}));
%             end
%         end
        
        for m = maxLev : 7;
            mipSPN = [mipDir num2str(m) '\'];
            mipTPN = [mipDir num2str(m+1) '\'];
            if ~exist(mipTPN,'dir'),mkdir(mipTPN);end
            dirMipSPN = dir([mipSPN '*.png']);
            inams = {dirMipSPN.name};
            
            clear r c
            for i  = 1:length(inams);
                nam = inams{i};
                A = sscanf(nam,'%d_%d');
                r(i) = A(1);
                c(i) = A(2);
            end
            
            colormap gray(255)
            
            parfor y = 0:ceil(max(r)/2)
                for x = 0:ceil(max(c)/2)
                    I = zeros(iSize,'uint8');
                    for yD = 0:1;
                        for xD = 0:1;
                            targ = find((r == (y * 2 + yD)) & (c ==  (x * 2 + xD)));
                            if ~isempty(targ)
                                Iraw = imread([mipSPN inams{targ}]);
                                Idown = imresize(Iraw,.5,'bilinear');
                                I(yD*iSize/2+1:(yD+1)*iSize/2, xD*iSize/2+1:(xD+1)*iSize/2) = Idown;
                            end
                        end
                    end
%                     image(I)
%                     pause(.1)
                    
                    newName = sprintf('%d_%d.png',y,x);
                    imwrite(I,[mipTPN newName]);
                end
            end
            
        end
        
        
    end
end