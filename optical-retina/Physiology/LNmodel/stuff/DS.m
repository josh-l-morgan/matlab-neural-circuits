%this creates plot of the temporal structure of the receptive field center
%first-pixel anomaly correction
T = [-500:10];
ROI = find (std(STA,0,2) >= 2*mean(std(STA,0,2)));
nROI = length(ROI);

mg = mean(STA(:)); low = min(STA(:)); high = max(STA(:));
STA(1,:) = mg;

hFig1 = figure('Position',[50 100 500 250]);
if (high-mg) > (mg-low)        %FOR ON CELLS
    ON = STA(ROI,:)-mg;
    [Y,I] = max(ON,[],2);
    ISort = sort(I);
    Early = mean(diag(I <= ISort(round(nROI/4))) * ON);
    Late = mean(diag(I >= ISort(round(nROI/4))...
        & I <= ISort(round(nROI/2))) * ON);
    plot(T, Early/max(Early), 'b', 'LineWidth',2)
    hold on
    plot(T, Late/max(Late), 'r', 'LineWidth',2)

else
    OFF = STA(ROI,:)-mg;    %FOR OFF ALIGNMENT
    [Y,I] = min(OFF,[],2);
    ISort = sort(I);
    Early = mean(diag(I <= ISort(round(nROI/4))) * OFF);
    Late = mean(diag(I >= ISort(round(nROI/4))...
        & I <= ISort(round(nROI/2))) * OFF);
    plot(T, -Early/min(Early), 'b', 'LineWidth',2)
    hold on
    plot(T, -Late/min(Late), 'r','LineWidth',2)

    ROIEarly = ROI(I <= ISort(round(nROI/4)));
    ROILate = ROI(I>= ISort(round(nROI/4))...
        & I <= ISort(round(nROI/2)));

end
axis tight
% figure(2)
% indicator1=zeros(4800,1);
% indicator1(median(ROIEarly)) = 0.5;
% indicator1(median(ROILate))=1;
% map1 = (reshape(indicator1,80,60))';
% pcolor(map1)
% shading flat
% colormap gray
% axis([22 58 12 48])

