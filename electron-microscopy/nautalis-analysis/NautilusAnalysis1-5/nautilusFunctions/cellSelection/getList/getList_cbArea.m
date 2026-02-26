function[colorList colPropRaw] = getList_cbArea(MPN);
    

   load([MPN 'cb2d.mat']);
    colorList = cb2d.IDs;
    colPropRaw = cb2d.areaUM;
