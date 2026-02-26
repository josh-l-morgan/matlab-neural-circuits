

MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export_14+04+25_matOut\'

load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


%%

I_c1001_Dim1 = showCellParts(obI,dsObj,1001,[1,1,1],1,fsize)*30;
image(I_c1001_Dim1)


I_c1001_Dim2 = showCellParts(obI,dsObj,1001,[1,1,1],2,fsize)*30;
image(I_c1001_Dim2)



I_c1001_Dim3 = showCellParts(obI,dsObj,1001,[1,1,1],3,fsize)*30;
image(I_c1001_Dim3)

I_1035 = showCellSum(obI,dsObj,1035,[1 0 0],Dim,fsize)*30;


I_1014 = showCellSum(obI,dsObj,1014,[0 1 0],Dim,fsize)*30;
