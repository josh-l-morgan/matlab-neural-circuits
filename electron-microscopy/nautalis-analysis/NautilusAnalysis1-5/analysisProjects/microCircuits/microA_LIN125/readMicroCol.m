
TPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\export\export_16+05+10_microALin125\';

textDir = dir([TPN '*.txt']);
    TFN = textDir(1).name;
    
fileName = [TPN TFN];


obI.colStruc = readVastColors(fileName);
obI.nameProps = getMicroNameProps(obI.colStruc.names);
