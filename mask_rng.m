function [xrng,yrng] = mask_rng(mask)
%MASK_RNG  Determines the range of a two-dimensional logical mask
%          within the matrix.
%
%          [XRNG,YRNG] = MASK_RNG(MASK) given the two-dimensional
%          logical mask, MASK, determines the range of the true values
%          within the matrix.  The minimum and maximum range in the X
%          (row) and Y (column) directions are returned in, XRNG and
%          YRNG, respectively.
%
%          NOTES:  None.
%
%          05-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in MASK_RNG:  An input is required!');
end
%
if ~islogical(mask)
  error(' *** ERROR in MASK_RNG:  Mask must be of type logical!');
end
%
npts = size(mask);
ndim = size(npts,2);
%
if any(npts<=1)||ndim~=2
  error([' *** ERROR in MASK_RNG:  Input logical mask must be', ...
         ' two-dimensional!']);
end
%
% Determine the Range of True Values
%
[xm,ym] = meshgrid(1:npts,1:npts);     % Generate mesh
xm(~mask) = NaN;       % Set area outside of mask to NaN
ym(~mask) = NaN;       % Set area outside of mask to NaN
%
xrng = [min(ym(:)) max(ym(:))];        % Image X range
yrng = [min(xm(:)) max(xm(:))];        % Image Y range
%
return