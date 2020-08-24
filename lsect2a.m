function [xycoords,idx] = lsect2a(m1,b1,pll)
%LSECT2A Finds all the intersections of a line with a piecewise linear
%        line in a 2-D plane.
%
%        XYCOORDS = LSECT2A(M1,B1,PLL) finds the intersections of
%        a line defined by a slope (M1) and intercept (B1) with a
%        piecewise linear line (PLL).  PLL is defined by a series of
%        2-D points with the X and Y coordinates of the points in
%        columns 1 and 2, respectively.  The X and Y coordinates of the
%        intersections are returned in the matrix XYCOORDS with the
%        X coordinates in the first column and the Y coordinates in the
%        second column.
%
%        [XYCOORDS,IDX] = LSECT2A(M1,B1,PLL) returns the index
%        IDX to the segment within PLL which had the intersections.
%
%        NOTES:  1.  Returns NaNs for coordinates if there is no
%                intersection.
%
%                2.  See lsect2.m which finds just the first
%                intersection.
%
%        20-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in LSECT2A: Three input arguments are required!');
end
%
% Check Vectors
%
[n1,l1] = size(m1);
[n2,l2] = size(b1);
[n3,l3] = size(pll);
%
if ((n1>1)||(l1>1)||(n2>1)||(l2>1)||(n3<2)||(l3~=2))
  error(' *** ERROR in LSECT2A:  Error in input arguments!')
end
%
% Loop through Piecewise Linear Line
%
ns = n3-1;              % Number of segments
%
xcoord = NaN(ns,1);
ycoord = NaN(ns,1);
idx = zeros(ns,1);
%
for k = 1:ns
   l = k+1;
   m2 = (pll(l,2)-pll(k,2))/(pll(l,1)-pll(k,1));
   b2 = pll(k,2)-m2*pll(k,1);
   xcoord(k) = (b2-b1)/(m1-m2);
   if ((xcoord(k)<pll(l,1))&&(xcoord(k)>pll(k,1)))|| ...
      ((xcoord(k)>pll(l,1))&&(xcoord(k)<pll(k,1)))
     ycoord(k) = m2*xcoord(k)+b2;
     idx(k) = k;
   else
     xcoord(k) = NaN;
   end
end
%
idnan = isnan(xcoord);
if all(idnan)
  xycoords = [NaN NaN];
  idx = [];
else
  idv = ~idnan;
  xycoords = [xcoord(idv) ycoord(idv)];
  idx = idx(idv);
end
%
return