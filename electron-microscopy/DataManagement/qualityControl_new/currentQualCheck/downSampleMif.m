%% Turn mosaics from waffer into viewable (centers and downsample) mosaics

%% Get Waffer folder information
TPN = GetMyDir;
if exist([TPN 'mif.mat'])
    load([TPN 'mif.mat'])
else
    mif = getMif(TPN);
end

if exist([TPN 'lowRes.mat'])
    load([TPN 'lowRes.mat'])
else
    lowRes.start = clock
end



for w = 1:length(mif.w)
    for s = 1:length(mif.w(w).sec)
        
        Info = imfinfo([mif.w(w).sec(s).tileFolders{1} '\' mif.w(w).sec(s).tileNams{1}]);
        
        pixelOverlap = 0;
        xs = Info.Width;
        ys = Info.Height;
        dSamp = xs/100;
        
        colormap gray(256)
        csize = length(dSamp:dSamp:xs-pixelOverlap);
        esize = length(1:dSamp:xs);
        
        idMos = mif.w(w).sec(s).idMos;
        
        sprintf('reading section %d of %d',s,length(mif.w(w).sec))
        mosDim = size(idMos);
        mos = zeros((mosDim-1)*csize+esize,'uint8');
        rowIds = mif.w(w).sec(s).rowIds;
        colIds = mif.w(w).sec(s).colIds;
        
        for t  = 1:length(mif.w(w).sec(s).tileNams)
            
            ystart = (rowIds(t)-1)*csize+1;
            xstart = (colIds(t)-1)*csize+1;
            if rowIds==mosDim(1),yLap = 0; ysize = esize; else yLap = pixelOverlap; ysize = csize; end
            if colIds==mosDim(2),xLap = 0; xsize = esize; else xLap = pixelOverlap; xsize = csize; end
            I = imread([mif.w(w).sec(s).tileFolders{t} '\' mif.w(w).sec(s).tileNams{t}],...
                'PixelRegion',{[1 dSamp ys-yLap],[1 dSamp xs-xLap]});
            mos(ystart:ystart+size(I,1)-1,xstart:xstart+size(I,2)-1) = 255-I;
            image(mos),pause(.1)
        end
        image(mos),pause(.1)
        lowRes.w(w).sec{s} = mos;
    end
    
end

save([TPN 'lowRes.mat'], 'lowRes')






