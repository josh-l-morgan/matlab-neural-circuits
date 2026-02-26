

SPN = 'Y:\Active\Data\AT\Morgan Lab\VAST_hxR\waf003_single_png\'
TPN = 'Y:\Active\Data\AT\Morgan Lab\VAST_hxR\waf003_single_png_single\'
if ~exist(TPN,'dir'), mkdir(TPN), end

% slice row col mip

dSPN = dir([SPN '*.png']);

%inams = cat(1,dSPN.name);



slices = 0:19;
cols = 2:3;
rows = 2:3;



filename = sprintf('%d_%d_%d_%d.png', slices(1),rows(1),cols(1),0);

I = imread([SPN filename]);
[ys xs zs] = size(I);

Ys = length(rows) * ys;
Xs = length(cols) * xs;
Zs = length(slices);

imDatAll = zeros(length(rows)*length(I(:,1)),length(cols)*length(I(1,:)),length(slices),'uint8');

for slice=slices
    disp(slice)
    rowit=0;
    imDatAll = zeros(length(rows)*length(I(:,1)),length(cols)*length(I(1,:)),'uint8');
    fileNameSection = sprintf('%d.tif',slice);
    for row=rows
        colit=0;
        for col=cols
            tempFileName=sprintf('%d_%d_%d_%d.png', slice,row,col,0);
            tempImDat=imread([SPN tempFileName]);
            xrange=rowit*length(I(:,1))+1:rowit*length(I(:,1))+length(I(:,1));
            yrange=colit*length(I(1,:))+1:colit*length(I(1,:))+length(I(1,:));
            imDatAll(xrange,yrange)=tempImDat;
            colit=colit+1;
        end
        rowit=rowit+1;
    end
    imwrite(imDatAll,[TPN fileNameSection])
end











