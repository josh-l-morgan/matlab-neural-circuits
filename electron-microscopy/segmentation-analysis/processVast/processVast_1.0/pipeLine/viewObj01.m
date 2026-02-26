function[oImage] = viewObj01(obj,downSamp);


if ~exist('downSamp','var')
    downSamp = 5;
end

oSubs = obj.subs;
minDim = [];
maxDim = [];
for i = 1:length(oSubs) %% clip off remaining zeros
    oSub = oSubs{i};
    if ~isempty(oSub)
        if isempty(minDim)
            minDim = min(oSub,[],1);
            maxDim = max(oSub,[],1);
        else
            minDim = min(minDim,min(oSub,[],1));
            maxDim = max(maxDim,max(oSub,[],1));
        end
    end
end

minDim = max([0 0 0],minDim-1);
maxDim = maxDim-minDim+1;

oImage.minDim = minDim;
oImage.maxDim = maxDim;

%% Create sum image

colormap gray(256)
xyDim = ceil([maxDim(1) maxDim(2)]/downSamp);
allSum = zeros([xyDim length(oSubs)]);

for i = 1:length(oSubs)
    disp(sprintf('Getting sums of object %d of %d',i,length(oSubs)))
    oSub = oSubs{i};
    
    if ~isempty(oSub)
        
        for m = 1:3
            oSub(:,m) = oSub(:,m)-minDim(m);
        end
        
        oSub = ceil(oSub/downSamp);       
        
        tempSum = zeros(xyDim);
        
        xyInd = sub2ind(size(tempSum),oSub(:,1),oSub(:,2));
        uniqueInd = unique(xyInd);
        if length(uniqueInd)>1
            histVal = hist(xyInd,uniqueInd);
            tempSum(uniqueInd) = histVal;
        else
            tempSum(uniqueInd) = 1;%length(xyInd);
        end
        allSum(:,:,i) = tempSum;
    end
end

maxSum = max(allSum,[],3).^.3;
image(maxSum*256/max(maxSum(:)))

oImage.sumObj = maxSum;

%% combine colors

ids = obj.ids;
typeLab = [];
colPal = hsv(1000);
for i = 1:length(ids)
    nam = ids{i,2};
    ID = ids{i,1};
    if ID>0;
%         if sum(regexp(lower(nam),'background'))
%             colTab(ID,:) = [0 0 0];
%         elseif sum(regexp(lower(nam),'rgc'))
%             colTab(ID,:) = [.1 .1 .1] + rand(1,3)*1;
%         elseif sum(regexp(lower(nam),'fil'))
%             colTab(ID,:) = [.1 .1 .1]  + rand(1,3)*1;
%         elseif sum(regexp(lower(nam),'unk'))
%             colTab(ID,:) = [.1 .1 .1]  + rand(1,3)*1;
%         elseif sum(regexp(lower(nam),'segment 52'))
%             colTab(ID,:) = [.1 .1 .1]+ rand(1,3)*1;
%         else
%             colTab(ID,:) = [.1 .1 .1]  + rand(1,3)*1;
%         end
        colTab(ID,:) = colPal(ceil(rand*size(colPal,1)),:);
    end
    
end

colTab(colTab>1) = 1;

oImage.colTab = colTab;

%% Paint Colors
col = zeros(size(allSum,1), size(allSum,2), 3);

isDat = squeeze(sum(sum(allSum,1),2));

for i = 1:size(allSum,3)
   
    if isDat(i);
         disp(sprintf('drawing %d of %d',i,size(allSum,3)))
        grabTemp = allSum(:,:,i);
          scaleUp = max(.04,min(1,200/max(grabTemp(:))));
        grabTemp = (grabTemp>0) * 100 + (grabTemp * scaleUp);
        parfor c = 1:3
            col(:,:,c) = (col(:,:,c) +  grabTemp * colTab(i,c));
        end
        pause(.01)
    end
    
end

image(uint8(col.^1))

oImage.colorOb = col;


