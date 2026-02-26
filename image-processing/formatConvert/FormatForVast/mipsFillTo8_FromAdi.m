

clear all
SPN = 'X:\Active\externalShare\dLGNs1\Josh_LGN_aligned_v3\Josh_LGN_aligned_v3\em\mip0\';
SPN = 'E:\LGNs1\AdiAlignment\em\mip0\';
TPN = 'E:\LGNs1\AdiAlign_MipMaps\'

%%Filter
fKern = ones(2,1);
fKern = fKern ./ sum(fKern(:));




secDir = dir(SPN);
secFold = {secDir([secDir.isdir]>0).name};

iSize = 2048;
blankI = zeros(iSize,'uint8');
for s = 1:length(secFold)
    nam = secFold{s};
    und = regexp(nam,'_');
    sec = str2num(nam(1:und-1));
    
    if ~isempty(sec)
        imageDir = [SPN secFold{s} '\'];
        inDir = dir([imageDir '*.png*']);
        iNams = {inDir.name};
        r = zeros(length(iNams),1);
        c = r;
        iEmpty = r;
        mipDir = sprintf('%s%d\\',TPN,sec);
        mipDir0 = [mipDir '0\'];
        if ~exist(mipDir0,'dir'),mkdir(mipDir0);end
        parfor i = 1:length(iNams);
            nam = iNams{i};
           
            ids  = sscanf(nam,'%d_w%d_Sec%d_Montage_tr%d-tc%d.png%s');
            iEmpty = ~isempty(regexp(nam,'_empty'));
            r = ids(4);
            c = ids(5);
             newName  = sprintf('%s%d\\0\\%d_%d.png',TPN,sec,r,c);
            if ~exist(newName,'file')
                
            if iEmpty
                I = blankI;
            else
                I = imread([imageDir nam],'png');
                I = imfilter(I,fKern);
            end
            
            %newName  = sprintf('%s%d\\0\\%d_%d.j2k',TPN,sec,r,c);
            newName  = sprintf('%s%d\\0\\%d_%d.png',TPN,sec,r,c);
            imwrite(I,newName);
            %imwrite(I,newName,'CompressionRatio',5)
            end
        end
        
                
        for m = 0 : 7;
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
                    newName = sprintf('%d_%d.png',y,x);
                    if ~exist(newName,'file')
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
                    
                    %newName = sprintf('%d_%d.png',y,x);
                    %newName = sprintf('%d_%d.j2k',y,x);
                    %imwrite(I,[mipTPN newName]);
                    %imwrite(I,[mipTPN newName],'CompressionRatio',5)
                    imwrite(I,[mipTPN newName])
                    end
                end
            end
            
        end
        
    end
end
