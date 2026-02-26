function[] = axAnaCB(TPN, DPN)

%Find cell body

colormap gray(255)



%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'


        load([TPN 'Dots.mat'])
    
        %%{
        d=dir([TPN 'mask']); %get number of files in directory
        d=d(3:size(d,1));

        clear I IM 
        IM=imread([TPN 'mask\' d(1).name]);
        [ys xs zs] = size(IM);
        IM=zeros(ys, xs, size(d,1),'uint8');
        for im=1:size(d,1)
            IM(:,:,im)=imread([TPN 'mask\' d(im).name]);
            PercentRead=im/size(d,1)*100
        end
        
        [y x z] = find3(IM==3);
        CBpos=[mean(y) mean(x) mean(z)];
        Dots.Im.CBpos=CBpos;
        %}
        
        CBpos=Dots.Im.CBpos; %open and scale
        CBpos(1:2)=CBpos(1:2)*.103; CBpos(3)=CBpos(3)*.3;
        Dpos=Dots.Pos; Dpos(:,1:2)=Dpos(:,1:2)*.103; Dpos(:,3)=Dpos(:,3)*.3;
        
        Dist2CB=dist(Dpos,CBpos);
        Dots.Dist2CB=Dist2CB;
        save([TPN 'Dots.mat'],'Dots')




