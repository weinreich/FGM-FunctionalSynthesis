%file meant to calculate Equation 3 from the grant proposal.  All but the
%acosd(), which can be imaginary! NOTE THAT THIS ULTIMATELY BECAME EQN 4 IN THE PAPER.
function [x] = Eqn3(Wx,Wy,Wxy);

x = log(Wxy/(Wx*Wy))/(-2*sqrt(log(Wx)*log(Wy)));
