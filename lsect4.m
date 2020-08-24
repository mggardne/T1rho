function [ips,ts,idc,ierr] = lsect4(p1,p2,pll,tol,clim)
%LSECT4   Finds the intersections of a closed line with a piecewise
%         linear line.
%
%         IPS = LSECT4(P1,P2,PLL) finds the intersections of a closed
%         line defined by an initial point (P1) and a final point (P2)
%         with a piecewise linear line (PLL).  PLL is defined by a
%         series of 2-D points with the X and Y coordinates of the
%         points in columns.  The X and Y coordinates of the
%         intersections with the piecewise linear line (PLL) are
%         returned in IPS.  IPS is empty if there is no intersection.
%
%         [IPS,TS,IDC] = LSECT4(P1,P2,PLL) returns the distance along
%         the lines to the intersection.  TS(1) is the normalized (0 to
%         1) distance along the first line to the intersection.  TS(2)
%         is the distance along the segment of the piecewise linear line
%         (PLL) with the intersection.  The index (IDC) is to the first
%         point of the segment in the piecewise linear line (PLL) that
%         intersects the input line.  If there is no intersection, TS
%         and IDC are empty arrays.
%
%         [IPS,TS,IDC,IERR] = LSECT4(P1,P2,PLL,TOL) sets IERR to true
%         if no intersection is found within tolerance (TOL).  Default
%         tolerance is 1e-8.
%
%         [IPS,TS,IDC,IERR] = LSECT4(P1,P2,PLL,TOL,CLIM) displays a
%         warning if the condition number of the matrix is greater than
%         a limit, CLIM.  This is usually due to parallel or nearly
%         parallel lines.  The default condition number limit is 1e+8.
%
%         NOTES:  1.  The lines are assumed not to be parallel.
%
%                 2.  The input points must be of length two (2).  The
%                 piecewise linear line must have two (2) columns.
%
%                 3.  The M-file lsect3.m must be in the current path
%                 or directory.
%
%                 4.  See lsect.m and lsect2.m for two-dimensional
%                 line (2-D) intersections.
%
%                 5.  See lisect.m and lisect2.m for three-dimensional
%                 line (3-D) intersections with infinite lines.  See
%                 lisect3.m for the intersection between two closed
%                 lines.
%
%         14-Nov-2016
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in  LSECT4:  LSECT4 requires three input', ...
         ' arguments.']);
end
%
if (nargin<4)||(isempty(tol))||(tol<=0)
  tol = 1e-8;
end
%
if (nargin<5)||(isempty(clim))||(clim<=1)
  clim = 1e+8;
end
%
% Check Points and Piecewise Linear Line (PLL)
%
p1 = p1(:);
p2 = p2(:);
[n1,l1] = size(p1);
[n2,l2] = size(p2);
[n3,l3] = size(pll);
%
if (n1~=2)||(l1~=1)||(n2~=2)||(l2~=1)||(n3<2)||(l3~=2)
  error(' *** ERROR in LSECT4:  Error in input arguments.')
end
%
% Initialize Arrays
%
nl = n3-1;              % Number of lines in piecewise linear line (PLL)
%
v1 = p2-p1;             % Direction vector (and length) for closed line
%
ips = NaN(nl,2);
ts = NaN(nl,2);
idc = NaN(nl,1);
ierr = true;
%
% Loop through Piecewise Linear Line
%
for k = 1:nl
%
   l = k+1;             % Next point
%
   p2 = pll(k,:)';                     % Point on line
   v2 = pll(l,:)'-p2;                  % Direction of line
%
% Check for Intersection
%
   [ip,t,ier] = lsect3(p1,v1,p2,v2,tol,clim);
   if ~isempty(ip)&&~ier
     ips(k,:) = ip';
     ts(k,:) = t';
     ierr = false;
     idc(k) = k;
   end
%
end
idv = ~isnan(idc);
idc = idc(idv);
ips = ips(idv,:);
ts = ts(idv,:);
%
return