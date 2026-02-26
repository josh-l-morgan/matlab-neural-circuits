function[se]=stder(Dat)
%%Calculates standard error

se=std(Dat(:))/length(Dat(:));