function[cellIDs properties] = getList_biconality();

%%Returns a list of cells near the seed cell and the number of manually
%%identified multitubular arrays in their cell body

cellProp = [ 106.0000    0.1649
  107.0000    0.0938
  108.0000    0.3934
  109.0000    0.5746
  110.0000    0.7746
  111.0000    0.1461
  112.0000    0.1898
  116.0000    0.4052
  117.0000    0.4131
  120.0000    0.6370
  123.0000    0.5003
  129.0000    0.7516
  130.0000    0.0598
  131.0000    0.2014
  134.0000    0.7818
  135.0000    0.1693
  137.0000    0.0852
  148.0000    0.3139
  156.0000    0.5813
  159.0000    0.7805
  162.0000    0.5060
  163.0000    0.5463
  169.0000    0.7876
  170.0000    0.4473
  201.0000    0.3385
  203.0000    0.7558
  205.0000    0.7415
  207.0000    0.4219
  210.0000    0.4090
  212.0000    0.5029
  213.0000    0.1409
  215.0000    0.6636
  216.0000    0.1255
  217.0000    0.3949
  218.0000    0.5438
  224.0000    0.4389
  232.0000    0.0645
  237.0000    0.4092
  267.0000    0.5688
  268.0000    0.5552
  273.0000    0.5744
  907.0000    0.6172]


useList = getList_tracedCells;

foundCell = zeros(size(cellProp,1),1)>0;
for i = 1:size(cellProp,1)
    foundCell(i) = sum(useList == cellProp(i,1))>0;
end
cellIDs = cellProp(foundCell,1);

properties = cellProp(foundCell,2);

