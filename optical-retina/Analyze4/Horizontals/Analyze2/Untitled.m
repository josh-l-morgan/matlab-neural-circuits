for k = 1 : size(Idir,1)
if ~exist([KPN 'Ter']), mkdir([KPN 'Ter']),end    
Names(k,1)={Idir(k).name};    
I=imread([TPN Idir(k).name]);
Ia(:,:,k)=I;

end
AbMax=max(Ia,[],3);
image(AbMax)



imwrite(AbMax,[KPN 'AbMax.tif'  ],'Compression','none')
