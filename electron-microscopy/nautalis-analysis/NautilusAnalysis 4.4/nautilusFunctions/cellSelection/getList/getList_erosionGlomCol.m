function[cellIDs cellProp] = getList_erosionGlomCol
colorTable = [ 0 1 0 ;  1 .2 0 ; 0 .2 1 ; 1 0 1 ; 1 1 0 ; 0 .2 1 ] * 1.5 ;


col = [];
allMembers = [];
useClade = [ 1 2 3 4 5];

for i = 1:length(useClade)
   members = getList_erosionGlom(useClade(i));
   
   %members = cladePick.members{useClade(i)};
   allMembers = [allMembers members];
   col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
    
end

cellIDs = allMembers;
cellProp = col;