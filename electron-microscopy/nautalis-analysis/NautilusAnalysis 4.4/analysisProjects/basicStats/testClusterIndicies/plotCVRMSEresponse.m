
%%
colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = zeros(1,n);
        biasedDist(1) = v;
        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,1)

image(cmat * 256/max(cmat(:)))

subplot(2,2,2)

image(rmat * 256/max(rmat(:)));


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = ones(1,n);
        biasedDist(1) = v;
        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,3)

image(cmat * 256/max(cmat(:)))

subplot(2,2,4)

image(rmat * 256/max(rmat(:)));

%%
clf


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = zeros(1,n);
        biasedDist = biasedDist + rand(1,n) * v;
        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,1)

image(cmat * 256/max(cmat(:)))

subplot(2,2,2)

image(rmat * 256/max(rmat(:)));


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = ones(1,n);
        biasedDist = biasedDist + rand(1,n)* v;
        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,3)

image(cmat * 256/max(cmat(:)))

subplot(2,2,4)

image(rmat * 256/max(rmat(:)));



%%
clf


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = zeros(1,n);
        biasedDist = biasedDist + (rand(1,n)*2).^v;
                biasedDist = biasedDist/sum(biasedDist);

        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,1)

image(cmat * 256/max(cmat(:)))

subplot(2,2,2)

image(rmat * 256/max(rmat(:)));


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = ones(1,n);
        biasedDist = biasedDist + (rand(1,n)*2).^v;
        biasedDist = biasedDist/sum(biasedDist);
        [c r] = cvrmse(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,3)

image(cmat * 256/max(cmat(:)))

subplot(2,2,4)

image(rmat * 256/max(rmat(:)));



%%
clf


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = zeros(1,n);
        biasedDist = biasedDist + (rand(1,n)*2).^v;
                biasedDist = biasedDist/sum(biasedDist);

        c = std(biasedDist)/mean(biasedDist);
        r = std(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,1)

image(cmat * 256/max(cmat(:)))

subplot(2,2,2)

image(rmat * 256/max(rmat(:)));


colormap jet(256)
clear cmat rmat

for n = 1:100
    for v = 1:100
        
        biasedDist = ones(1,n);
        biasedDist = biasedDist + (rand(1,n)*2).^v;
        biasedDist = biasedDist/sum(biasedDist);
        %[c r] = cvrmse(biasedDist);
        c = std(biasedDist)/mean(biasedDist);
        r = std(biasedDist);
        cmat(n,v) = c;
        rmat(n,v) = r;
        
    end
end

subplot(2,2,3)

image(cmat * 256/max(cmat(:)))

subplot(2,2,4)

image(rmat * 256/max(rmat(:)));




