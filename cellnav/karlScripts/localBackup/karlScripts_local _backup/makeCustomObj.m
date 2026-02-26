function objList=makeCustomObj(includeStrings,excludeStrings,tis)
% This will make you a custom object (maybe also FV) from segmentation
% objects with the inclusion and exclusion strings specified.
% Example:
% makeCustomObj({["1216","syn"],["1197","syn"]},{["probably"]},tis)
% ( ( 1216 & syn ) | ( 1197 & syn ) ) & ~probably
% includes '1216' AND 'syn', or '1197' AND 'syn', but nothing with 'probably'
% good luck

