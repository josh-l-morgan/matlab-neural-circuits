%function[focSec] = checkFocus(wif)

%% Turn mosaics from waffer into viewable (centers and downsample) mosaics
clear all 

%% Define variables
dsampy = 100; %down sample image reading
dsampx = 100;
        
yshift = [ 0 0 0 1 1 1 2 2 2]; %shifts for contrast measurement
xshift = [ 0 1 2 0 1 2 0 1 2];



%% Get Waffer folder information
wif = GetMyWafer;
%% Parse XML
[tree, rootname, dom]=xml_read(wif.xml{end});

xs = tree.MosaicSetup.TileWidth;
ys = tree.MosaicSetup.TileHeight;

%% Target Dir
TPN = wif.dir; TPN = [TPN(1:end-1) 'Shaped\'];
TPNsav = [TPN 'quality\'];
if ~exist(TPNsav),mkdir(TPNsav);end

%% test Focus
colormap gray(256)


for s = 1 : length(wif.sec) % run sections
    sprintf('reading  section %d of %d',s,length(wif.sec))

    %[tree, rootname, dom]=xml_read(wif.sec(s).xml);
     rc = wif.sec(s).rc;
     mosDim = max(rc,[],1);
     mos = zeros(mosDim(1),mosDim(2),3);
    for t = 1:length(wif.sec(s).tile)
        
        I1 = double(imread(wif.sec(s).tile{t},'PixelRegion',{[ 1 dsampy ys],[1 dsampx xs]}));
        %I1 = I1(1:end-1,:);
        Id = zeros(size(I1,1),size(I1,2),length(yshift));
        
        for i = 1 :length(yshift)
            Is = double(imread(wif.sec(s).tile{t},'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
            Is = Is(1:size(I1,1),1:size(I1,2));
            Id(:,:,i) = abs(Is-I1);
            trackId(i,1:2) = [yshift(i) xshift(i)];
        end
        %Ic(:,:) = min(abs(Is-circshift(I1,[5 0])),abs(Is-circshift(I1,[0 5])));
        Ic(:,:) = abs(Is - circshift(I1,[23 0]));
        vals = sort(Ic(:),'descend');
        difC = median(vals(:));

%%     Find contrasts
        % % %
        % % %
        % % %
          
        I = double(imread(wif.sec(s).tile{t},'PixelRegion',...
            {[ round(ys/2) - dsampy round(ys/2)+dsampy],[round(xs/2)- dsampx round(xs/2)+dsampx]}));
        subplot(1,2,1)
        image(I)
        
        fronts = {[1 4 7 2 5 8],[1 2 3 4 5 6]};
        backs = {[2 5 8 3 6 9],[4 5 6 7 8 9]};
        comp1 = {[1 2 3],[1 2 3]};
        comp2 = {[4 5 6],[4 5 6]};
         
        for f = 1: length(fronts)
            difI = mean(Id(:,:,fronts{f})) - mean(Id(:,:,backs{f}));
            ddI = difI(:,:,comp1{f}) - difI(:,:,comp2{f});
            mdI(:,:,f) = abs(mean(ddI,3));
        end
        
        %%Center
        Icenter = Id(:,:,5) - (sum(Id,3) - Id(:,:,5))/8;
        mIcenter = mean(Icenter(:))
        %% xdim
        vals = mdI(:,:,1);
        sIx = vals(:)/mean(vals(:));
        sIx = sort(sIx,'descend');
        subplot(2,2,2)
        hist(sIx)
        xlim([0 12])
        ylim([0 200])
        
        %% ydim
        vals = mdI(:,:,2);
        sIy = vals(:)/mean(vals(:));
        sIy = sort(sIy,'descend');
        subplot(2,2,4)
        hist(sIy)
        
        xlim([0 12])
        ylim([0 200])
        pause
        
        difx = mean(sIx(1:fix(length(sIx)/300)));
        dify = mean(sIy(1:fix(length(sIx)/300)));
        con = difC;
        
        dmap = [difx dify];
        mos(rc(t,1),rc(t,2),:) = [dmap(1) dmap(2) con];
        
%% Make contrast image
        freqMap = sum(Id,3);
        bKern = ones(4,10);
        freqReg = fastCon(freqMap,bKern);
        [h w] = find(freqReg == max(freqReg(:)),1);
        focTarg = [h/size(freqMap,1) w/size(freqMap,2)];
        
        focSec(s).tile(t).difX1X2Y1Y2XY1 = [mean(sIx) mean(sIy)];
        focSec(s).tile(t).focY = dify;
        focSec(s).tile(t).focX = difx;
        focSec(s).tile(t).contrast = con;
        focSec(s).tile(t).focusTarget = focTarg;
%%
        
    end
    mosA(:,:,:,s) = mos;
    
end
save([TPNsav 'qual.mat'],'mosA')
save([wif.dir 'focSec.mat'],'focSec')









