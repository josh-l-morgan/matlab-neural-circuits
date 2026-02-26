function[CB]=mosaic(Num);
%Num=100
%%Distributes Num number of nodes evenly across a 2 d field 
%%With the resolution of fsize

colormap gray(256)
Grad=1;
CB.pos=rand(Num,2);

fsize = 200
reps = 50
CB.pos=CB.pos*fsize;
CB.pos=round(CB.pos);
CB.pos(CB.pos<2)=2;
CB.pos(CB.pos>(fsize-2))=fsize-2

%Define filter
    G1=fspecial('gaussian',fsize,fsize/(20));
    G1=G1* 1/max(G1(:));
    G2=fspecial('gaussian',fsize,fsize/(5));
    G2=G2* 1/max(G2(:));
    
    G=G1;
    


for rep = 1:reps
    pdone=rep/reps*50
    %scoot from wall
    CB.pos(CB.pos<2)=2;
    CB.pos(CB.pos>(fsize-2))=fsize-2;

    field = zeros(fsize,fsize);

    for i = 1: Num
        field(CB.pos(i,1),CB.pos(i,2))=10;
    end




    R=imfilter(field,G);
    
    
    
%        
%     Wall([1 fsize],1:fsize)=10;
%     Wall(1:fsize,[1 fsize])=10;
%     Gw=fspecial('gaussian',fsize,fsize/10);
%     Gw=Gw * 1/ max(Gw(:));
%     Rw=imfilter(Wall,Gw);
%     
%     R=R+Rw*Grad;
    
    %image(R*255/max(R(:))) ,pause(.1)
   

    for n = 1: Num

        Near=R(CB.pos(n,1)-1:CB.pos(n,1)+1,CB.pos(n,2)-1:CB.pos(n,2)+1);
        [miny minx]=find(Near==min(Near(:)));
        
        if min(Near(:))/R(CB.pos(n,1),CB.pos(n,2))<.99
        Targ=fix(rand*size(miny,1))+1;
        CB.pos(n,1)=CB.pos(n,1)+miny(Targ)-2;
        CB.pos(n,2)=CB.pos(n,2)+minx(Targ)-2;
        end

    end %move nodes
     image(field*255/max(field(:))) ,pause(.01)
    
end %reps


%% Show CBs
field = zeros(fsize,fsize);

for i = 1: Num
    field(CB.pos(i,1),CB.pos(i,2))=10;
end

R=imfilter(field,G);

image(R*255/max(R(:))) ,pause(.1)
