function [ absresdiff ] = fitabsres3( k3, doseavg, res3avg )
%This function returns a function to be minimized by lsqnonlin. The
%function to be returned is the differnece between the model fit for the
%absolute value of the residuals and the experimental data of the two
%population model

slopelin = k3(1);
expterm = k3(2);
shift = k3(3);

%y = 5*x.*exp(-.05*x);

absresdiff = (slopelin*(doseavg+shift).*exp(-expterm*doseavg)) - res3avg;


end