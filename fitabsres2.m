function [ absresdiff ] = fitabsres2( k2, doseavg, res2avg )
%This function returns a function to be minimized by lsqnonlin. The
%function to be returned is the differnece between the model fit for the
%absolute value of the residuals and the experimental data of the two
%population model

slopelin = k2(1);
expterm = k2(2);
shift = k2(3);

%y = 5*x.*exp(-.05*x);

absresdiff = (slopelin*(doseavg+shift).*exp(-expterm*doseavg)) - res2avg;


end
