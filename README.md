# matrix-information-metrics
Mutual information and matrix MMSE, for a matrix SNR input

INFO_function(SList,X,p)
% Computes the mutual information I_X(S) and the MMSE matrix  M_X(s) for a
% discrete random vector supported on the columns of X with pmf p, SNR is a variable.

function [I,M] = INFO_function_2GMM(s,aXList,mucell,sz)
% This function compute the mutual information I_X(s) and the MMSE M_X(s)
% for a Gaussian mixture with paramters aXList, muXList, and tauX2List. s is the matrix SNR.

function [I,M] = INFO_function_GMM(sList,aXList,muXList,tauX2List)
% This function compute the mutual information I_X(s) and the MMSE M_X(s)
% for a Gaussian mixture with paramters aXList, muXList, and tauX2List.  SNR is a variable.

Entropy_3comm_corrected(aXList,mucell)
Compute the entropy for a GMM with 3 clusters

