


SPN = 'Z:\Data\hxH\export\';


inams = dir([SPN '*.png']);
sub = [];
for i = 1:length(inams)
    disp(sprintf('reading %d of %d',i,length(inams)))
    I = imread([SPN inams(i).name]);
    [y x] = find(I>0);
    v = I(I>0);
    z = ones(length(y),1) * i;
    sub = cat(1,sub,[y x z]);
end


renderOb = 1;
col = [1 0 0];

    downSamp = 8;
    if exist('crop','var')
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
        sub = sub(useSub,:);
        
    end
    smallSub = shrinkSub(sub,downSamp);
    tic
    if isempty(smallSub)
        disp(sprintf('no points on %d',cellList{i}))
    else
    fv = subVolFV(smallSub,[],renderOb);
    renderFV(fv,col);
    pause(.01)
    hold on
    fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
    fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
    %STLWRITE(FILE, FACES, VERTICES)
    %stlwrite(fileNameSTL,fv.faces,fv.vertices);
    vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
    toc
    %     cellDat(i).subs = sub;
    %     cellDat(i).fv = fv;
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList{i},i,length(cellList)));
