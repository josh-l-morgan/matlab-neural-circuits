%%Create mosaics using simple distance as probability

fsize=500;
field=zeros(fsize);
Num=20;
gx=rand(Num,1)*fsize;
gy=rand(Num,1)*fsize;
Dmap=zeros(fsize,fsize,Num);
Dfac=1;  %exponent to which distance is raised to calculate probability 1, or 2 is normal. 

for y = 1: fsize
    for x = 1:fsize

        Dist=sqrt((gy-y).^2+(gx-x).^2);
        Dmap(y,x,:)=Dist;
        Prob=1./((Dist/fsize).^Dfac);
        Prob=Prob/sum(Prob);
        Pick=rand;
        psum=0;
        for i = 1:Num
            psum=psum+Prob(i);
            if psum>Pick
                Targ=i; %pick cell
                break
            else
                Targ=Num;
            end
        end
        field(y,x)=Targ; %color field
        
    end
    PercentDoneDrawing=y/fsize *100
end

colormap hsv(256)
subplot(2,1,1)
image(field*(256/Num)),pause(.01)

%% Bin 
Dens=zeros(Num,fsize/2);
for i = 1:Num  
    for r =1:fsize/2;
       Mask=Dmap(:,:,i)>(r-5) & Dmap(:,:,i)<(r+5);
       Self(r)=sum(field(Mask)==i);
       All(r)=sum(Mask(:));
    end
    Dens(i,:)=Self./All;
    
end


%% Normalize Dens
subplot(2,1,2)
for i = 1:Num
    dNorm(i)=min(find(Dens(i,:)<.2));
    nDens=interp(Dens(i,:),10);
    nDens=imresize(nDens,[1,1000/dNorm(i)]);
    plot(nDens,'r'); pause(.01)
    hold on
end
hold off













