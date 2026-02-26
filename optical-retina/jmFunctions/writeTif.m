function[] = writeTif(I,name);

TPN = GetMyDir;
imwrite(I,[TPN name],'Compression','none');