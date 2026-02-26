

fileName = 'D:\LGNs1\Segmentation\VAST\S8\joshm\S8_joshmHome_14+04+17.vss.txt';

obI.colStruc = readVastColors(fileName);
obI.nameProps = getNameProps(colStruc.names);

cellIDs = obI.colStruc.cell;