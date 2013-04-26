%file meant to calculate fitness
function [fitness] = fitness(Selection,x);

fitness = exp(-0.5*x'*Selection*x);
