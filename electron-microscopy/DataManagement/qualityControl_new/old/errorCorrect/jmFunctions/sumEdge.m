function[sI] = sumEdge(I,threshs, sizes,gausRad);
sI = I * 0;

if ~exist('threshs')
    threshs= [.3:.1:.7]
end

if ~exist('sizes')
    sizes = 0.5:0.5:5
end


for t = threshs
    for s = sizes
    thresh = t/ (s/min(sizes));
    sI = sI + edge(I,'Canny',[.001 thresh],s);
    image(fitH(sI)),pause(.01)
    %[t s]
%     pause
    end
end

%if exist('gausRad')
if ~exist('gausRad')
    gausRad = 1;

    gKern = gaus3d([gausRad * 3 gausRad * 3 1],gausRad);
    sI = fastCon(sI,gKern);
elseif gausRad
    gKern = gaus3d([gausRad * 3 gausRad * 3 1],gausRad);
    sI = fastCon(sI,gKern);
end
image(sI* 50/mean(sI(sI>0)))

%end
