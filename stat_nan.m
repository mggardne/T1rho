function [rmean,sd,n] = stat_nan(data,idnn,idv)
%STAT_NAN  Given a data matrix, a linear index matrix to the nine (9)
%          nearest neighbors and a logical matrix to the valid indices
%          in the index matrix, returns the means and optional standard
%          deviations of the data matrix as indexed by the columns of
%          the linear index matrix.
%
%          RMEAN = STAT_NAN(DATA,IDNN,IDV) given a data matrix, DATA,
%          a linear index matrix to the nine (9) nearest neighbors,
%          IDNN, and a logical matrix to the valid indices in the index
%          matrix, IDV, returns the means, RMEAN, of the data matrix,
%          DATA, as indexed by the columns of the linear index matrix,
%          IDNN.
%
%          [RMEAN,SD] = STAT_NAN(DATA,IDNN,IDV) returns the sample
%          standard deviations, SD, of the data matrix, DATA, as
%          indexed by the columns of the linear index matrix, IDNN.
%
%          [RMEAN,SD,N] = STAT_NAN(DATA,IDNN,IDV) returns the number of
%          observations, N, of the data matrix, DATA, as indexed by the
%          columns of the linear index matrix, IDNN.
%
%          NOTES:  1.  For two-dimensional (2D) matrices.
%
%                  2.  If a column of the linear index matrix does not
%                  have at least two valid indices, the mean or
%                  standard deviations will return as Inf or NaN.
%
%                  3.  For use with get_nn_idx.m to calculate the means
%                  and standard deviations of the nearest neighbors.
%                  See get_nn_idx.m.
%
%          11-Sep-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  error(' *** ERROR in STAT_NAN:  Three inputs are required!');
end
%
% Convert NaNs to 1 and Use Linear Index and Logical Matrices to
% Calculate the Means
%
idnn(~idv) = 1;
%
dat = data(idnn);
dat = dat.*idv;         % Set NaN indices to zeros
n = sum(idv);           % N
rmean = sum(dat)./n;
%
% Calculate the Standard Deviations
%
if nargout>1
  rmat = repmat(rmean,size(idnn,1),1); % Matrix of means
  dat = dat-rmat;
  dat = dat.*idv;       % Set NaN indices to zeros
  dat = sum(dat.*dat);  % Sum of squares
  sd = sqrt(dat./(n-1));               % Sample SD
end
%
return