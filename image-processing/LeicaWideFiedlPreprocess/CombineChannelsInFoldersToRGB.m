

SPN = 'G:\RReSTOReMeeting\';
SPN = 'G:\RReSTOReMeeting\SM25-1\';
SPN = 'G:\RReSTOReMeeting\SMjlm2_26-1\';
SPN = 'G:\RReSTOReMeeting\LxN_Retina_Cage9M1_\';
SPN = 'G:\RReSTOReMeeting\LxN_Retina_Cage9M0_PTEN+CNTF\';
SPN = 'Z:\Active\morganLab\DATA\LGN_Regeneration\LxP\';
SPN = 'Z:\Active\morganLab\DATA\LGN_Regeneration\rrestoreMeeting\rrestoreLGNs\LxL_CNTF_c3c4\';
SPN = 'Z:\Active\morganLab\DATA\LGN_Regeneration\rrestoreMeeting\rrestoreother\LxN_OpticTracts\';
SPN = 'Z:\Active\morganLab\DATA\LGN_Regeneration\rrestoreMeeting\rrestoreother\LxN_OpticTractsMore\';
SPN = 'H:\';

dSPN = dir([SPN '*.tif']);
iNams = {dSPN.name}';
iNum = length(iNams);
useChan = [3 2 4];
chanCon = [1 1 1];
chanOffset = [-5 -5 -5];


clear base chan
for i = 1:length(iNams)
    nam = iNams{i};
    raw = regexp(nam,'_RAW');
    if isempty(raw)
        base{i,1} = [];
    else
        base{i,1} = nam(1:raw(1)-1);
    end

    ch = regexp(nam,'_ch');
    if isempty(ch)
        chan(i,1) = 0;
    else
        channelNumber = str2num(nam(ch(end)+3:ch(end)+4));
        if isempty(channelNumber)
            chan(i,1) = 0;
        else
            chan(i,1) = channelNumber;
        end
    end
end

baseID = zeros(iNum,1);
for i = 1:length(iNams)
    for o = 1:length(iNams)
        if strcmp(base{i},base{o})
            baseID(i,1) = o;
        end
    end
end

baseIDs = setdiff(unique(baseID),0);

for i = 1:length(baseIDs)

    isBase = baseID==baseIDs(i);
    chID(1) = find(isBase & (chan==useChan(1)));
    chID(2) = find(isBase & (chan==useChan(2)));
    chID(3) = find(isBase & (chan==useChan(3)));

    clear I
    for c = 1:3
        Ic = double(imread([SPN iNams{chID(c)}]));
        
        if 0
            
            Itweak = Ic * chanCon(c);
            Itweak = Itweak * 256 / (2^16);
            Itweak = Itweak + chanOffset(c);
        else
            Im = medfilt2(Ic,[3 3]);
            maxMed = max(Im(:));
            Itweak = Ic * 256/maxMed;
        end
        I(:,:,c) = Itweak;
    end

    fileName = sprintf('%s%s_RGBmax.tif',SPN,base{baseIDs(i)});
    imwrite(uint8(I),fileName);
    image(uint8(I))
    pause(1)

end





