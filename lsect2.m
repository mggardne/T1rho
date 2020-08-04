function [xcoord,ycoord,idx] = lsect2(m1,b1,pll)
%LSECT2  Finds the first intersection of a line with a piecewise linear
%        line in a 2-D plane.
%        [XCOORD,YCOORD] = LSECT2(M1,B1,PLL) finds the intersection of
%        a line defined by a slope (M1) and intercept (B1) with a
%        piecewise linear line (PLL).  PLL is defined by a series of
%        2-D points with the X and Y coordinates of the points in
%        columns 1 and 2, respectively.
%
%        [XCOORD,YCOORD,IDX] = LSECT2(M1,B1,PLL) returns the index IDX
%        to the segment within PLL which had the intersection.
%
%        NOTES:  1.  Returns NaNs for coordinates if there is no
%                intersection.
%
%        25-Sep-96 * Mack Gardner-Morse
%
%        25-Jun-2009 * Mack Gardner-Morse * Updated to return index IDX.
%
%        06-Apr-2020 * Mack Gardner-Morse * IDX is empty if no
%                                           intersection.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error('LSECT requires three input arguments.');
end
%
% Check Vectors
%
[n1,l1] = size(m1);
[n2,l2] = size(b1);
[n3,l3] = size(pll);
%
if ((n1>1)||(l1>1)||(n2>1)||(l2>1)||(n3<2)||(l3~=2))
  error('LSECT2:  Error in input arguments.')
end
%
% Loop through Piecewise Linear Line
%
ycoord = NaN;
%
for k = 1:n3-1
   l = k+1;
   m2 = (pll(l,2)-pll(k,2))/(pll(l,1)-pll(k,1));
   b2 = pll(k,2)-m2*pll(k,1);
   xcoord = (b2-b1)/(m1-m2);
   if ((xcoord<pll(l,1))&&(xcoord>pll(k,1)))||((xcoord>pll(l,1))&&(xcoord<pll(k,1)))
     ycoord = m2*xcoord+b2;
     break;
   end
   xcoord = NaN;
end
%
if isnan(xcoord)
  idx = [];
else
  idx = k;
end
%
return