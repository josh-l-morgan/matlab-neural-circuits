clear all
colormap gray(256)
N = 1000000;
fsize = 200;
edge = 0.95;
check1 = 0.96;
check2 = 0.97;
difRate = .1;

vol = zeros(fsize,fsize,fsize);
p = rand(N,3)* fsize;
p(p>fsize) = 1;
p(p<1) = fsize;
p(:,1) = p(:,1)* edge;

for i = 1:10000
    
    m = (rand(N,3)-.5) * difRate;
    p = p + m;
    p(p>fsize) = 1;
    p(p<1) = fsize;


    rec1(i) = sum(round(p(:,1))==round(fsize*check1));
    rec2(i) = sum(round(p(:,1)) == round(fsize*check2));

    if ~mod(i,100)
        
        rp = fix(p) + 1;
        rp(rp>fsize) = 1;
        iVol = sub2ind([fsize fsize fsize],rp(:,1),rp(:,2),rp(:,3));
        vol = vol * 0;
        vol(iVol) = 1;
        show = 255 - (sum(vol,3) * 100);

        subplot(1,2,1)
        image(fitH(show)); pause(.01)

        subplot(1,2,2)
        plot(rec1,'b')
        hold on
        plot(rec2,'r')
        hold off
    end


end
