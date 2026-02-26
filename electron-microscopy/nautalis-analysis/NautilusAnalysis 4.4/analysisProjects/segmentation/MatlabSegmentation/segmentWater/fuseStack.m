clear all
TPN = GetMyDir;

% Inams = getPics([TPN]);
% if exist([TPN 'edited'])
%     labelDir = 'edited\';
% elseif exist([TPN 'labeled'])
%     labelDir = 'labeled\';
% else
%     sprintf('No labels found.')
% end

labelDir = 'labeled\';
inams = getPics([TPN labelDir]);  %find all tifs

if ~exist([TPN 'fused'])
    mkdir([TPN 'fused']);
end

I1 = imread([TPN labelDir inams{1}]);
rProps =regionprops(I1,'Area');
Areas(:,1) = [rProps.Area]';
Areas(:,2) = Areas(:,1) * 0 + 1;

myCol = hsv(256);
myCol(:,1) = 0;
colormap colorcube(256)
image(mod(I1,256)),pause(.01)

id16 = 1:65535; %list usable ids
id16 = setdiff(id16,I1(:));

mergeThresh = .4;
for i = 2:length(inams)

    sprintf('running plane %d of %d',i,length(inams))
    maxI1 = max(I1(:));
    sprintf('Ids are at %d',maxI1)
    I2 = imread([TPN labelDir inams{i}]);
    if length(unique(I2))<4  %crappy failsafe, replace with mask
        I2 = I1;
    end
    I2(I2>0) = id16(I2(I2>0));

    for r = 1:1
        I2props = regionprops(I2,'Area');
        iAreas = [I2props.Area];
        allIds = find(iAreas>0);
        iAreas = iAreas(iAreas>0);

        for a = 1: length(allIds)
            %sprintf('running area %d of %d',a,length(allIds))
            overlapped = I1(I2 == allIds(a));
            overlapped = overlapped(overlapped>0);
            if ~isempty(overlapped)
                [mO F] = mode(double(overlapped)); %apply threshold
                if F>= (iAreas(a) * mergeThresh)
                    I2(I2==allIds(a)) = mO;
                end
            end

        end  % end run merge forward

        %% run murge backward
        I1props = regionprops(I1,'Area');
        iAreas = [I1props.Area];
        allIds = find(iAreas>0);
        iAreas = iAreas(iAreas>0);

        for a = 1: length(allIds)
            %sprintf('running area %d of %d',a,length(allIds))
            overlapped = I2(I1 == allIds(a));
            overlapped = overlapped(overlapped>0);
            if ~isempty(overlapped)
                [mO F] = mode(double(overlapped)); %apply threshold
                if F>= (iAreas(a) * mergeThresh)
                    I1(I1==allIds(a)) = mO;
                    I2(I2 == allIds(a)) = mO;
                end
            end
        end

%% Display
        subplot(1,2,1)
        image(mod(I1,256))
        subplot(1,2,2)
        image(mod(I2,256))
        pause(.01)
    end
if length(id16)<10000
    id16 = 1:65500;
end
    id16 = setdiff(id16,[I2(:); I1(:)]);
    imwrite(uint16(I1),[TPN 'fused\fused' inams{i-1}],'Compression','none');
    I1 = I2;

end
imwrite(uint16(I2),[TPN 'fused\fused' inams{end}],'Compression','none');


%% Display

inams = getPics([TPN 'fused']);  %find all tifs
subplot(1,1,1)

cmap = colormap(hsv(256));
cmap(1,:) = 0;
colormap(cmap);

for i = 1: length(inams)
    I = imread([TPN 'fused\' inams{i}]);
    Imod = mod(I,256);
    Imod(~I) = 0;
    sprintf('showing plane %d',i)
    image(Imod)
    pause(1)
end

%%  PlotSkel
% %Make all Seg
% 'Display results'
% subplot(1,1,1)
% hold off
% plot(1,1)
% whitebg('k')
% 'render cell'
% uCol = hsv(100);
% skipdist = 11 + round(rand  * 8);
% for w = 2:max(storeI(:))
%     wInd = storeI == w;
%     if length(wInd) >10
%         subplot(1,1,1)
%         w
%         'rendering'
%     p = patch(isosurface(storeI==w,.1));
%     set(p,'FaceColor',uCol(mod(w*skipdist,100)+1,:),'EdgeColor','none');
%     daspect([1/4 1/4 1/30])
%     alpha(.8)
%     view(3); axis tight
%     camlight
%     lighting gouraud
%     'rendered'
%     pause
%     end
% end


%% Possible methods
%{
Thresholds
%}
