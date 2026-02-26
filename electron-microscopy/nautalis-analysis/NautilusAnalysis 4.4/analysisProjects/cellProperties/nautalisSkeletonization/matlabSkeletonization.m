


    testVol = zeros(max(surfVox.subs,[],1));
    testVol(sub2ind(size(testVol),surfVox.subs(:,1),surfVox.subs(:,2),surfVox.subs(:,3)))= 1;
    sk  = bwskel(logical(testVol));
    image(sum(sk,3)*100)