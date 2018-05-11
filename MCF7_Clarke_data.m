% Script takes literature values from 1991 Clarke Kaplan paper and plots
% the homogenous resistance to 2DG versus the homogenous resistance to Dox
% via the IC50 value.
clear all, close all, clc
IC50dox = [ 1 9 108 144 88]
IC502DG = [ 1 1.01 0.06 0.03 0.04]
plot ((IC50dox), (IC502DG), 'o')
xlabel ('IC50dox')
ylabel ('IC502DG')
title ('Plot of Dox Resistance vs. 2DG Toxicity')