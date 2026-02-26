function [err] = sem(input)

err = std(input)/sqrt(length(input));