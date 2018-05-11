function [ absresdiff ] = fitabsres1( k1, doseavg, res1avg )
%This function returns a function to be minimized by lsqnonlin. The
%function to be returned is the differnece between the model fit for the
%absolute value of the residuals and the experimental data of the bulk
%model

slopelin1 = k1(1);
expterm = k1(2);
shift = k1(3);

%y = 5*x.*exp(-.05*x);
%exp(slope1.*(dose - cen1)))
absresdiff = (slopelin1*(doseavg+shift).*exp(-expterm*doseavg)) - res1avg;


end
