TPN = GetMyDir

    TPNd=dir(TPN);TPNd=TPNd(3:length(TPNd));
    Names={};
    for i = 1:length(TPNd)
        nam=TPNd(i).name;
        ln=length(nam);
        if nam(ln-3:ln)=='.tif'
            Names{length(Names)+1}=nam;
        end
    end
    
    
    for i = 1: length(Names)
        I = imread([TPN Names{i}]);
        