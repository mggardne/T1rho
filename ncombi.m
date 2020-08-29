function [l,m] = ncombi(idx,nszc)
%NCOMBI  Finds the two indices into the combinatorial matrix NCOMB in
%        M-file T1rho_maps.m.  L is the number of elements in the
%        combinations minus one and M is the index into the combinations
%        with L+1 elements.
%
%        [L,M] = NCOMBI(IDX,NSZC) Finds the indices L and M into the
%        combinatorial matrix NCOMB given the index into the total
%        number of combinations IDX and the cumulative sum of the
%        number of combinations for each number of elements in the
%        combinations NSZC.
%
%        NOTES:  1.  See T1rho_maps.m for the construction of
%                combinatorial matrix NCOMB.
%
%        26-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in NCOMBI:  Two inputs are required!');
end
%
% Find Indices
%
l = find(idx>nszc,1,'last');
m = idx-nszc(l);
%
return