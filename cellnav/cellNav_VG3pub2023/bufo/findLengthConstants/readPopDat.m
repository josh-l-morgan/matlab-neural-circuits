



global glob

SPN =  [glob.datDir 'Analysis\Data\preproc\'];

datName = 'Circle_Polarity_Depth_MoveCircle072917.mat'
popDat = load([SPN datName])


scatter(popDat.Depth,popDat.Polarity)








