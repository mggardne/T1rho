function [ips,ts,idc,ierr] = lsect5(pl1,pl2,tol,clim)
%LSECT5   Finds intersections between two piecewise linear lines.
%
%         IPS = LSECT5(PL1,PL2) finds the intersections of a piecewise
%         linear line (PL1) with a second piecewise linear line (PL2).
%         PL1 and PL2 are defined by a series of 2-D points with the X
%         and Y coordinates of the points in columns.  The X and Y
%         coordinates of the intersections between the piecewise linear
%         lines are returned in IPS.  IPS is empty if there are no
%         intersections.
%
%         [IPS,TS,IDC] = LSECT5(PL1,PL2) returns the distance along
%         the lines to the intersection.  TS(1) is the normalized (0 to
%         1) distance along a segment in the first piecewise linear
%         line (PL1) to the intersection.  TS(2) is the distance along
%         the segment of the second piecewise linear line (PL2) to the
%         intersection.  IDC(1) is the index to the first point of the
%         segment in the first piecewise linear line (PL1) with an
%         intersection.  IDC(2) is the index to the first point of the
%         segment in the second piecewise linear line (PL1) with an
%         intersection.  If there is no intersection, TS and IDC are
%         empty arrays.
%
%         [IPS,TS,IDC,IERR] = LSECT5(PL1,PL2,TOL) sets IERR to true
%         if no intersection is found within tolerance (TOL).  Default
%         tolerance is 1e-8.
%
%         [IPS,TS,IDC,IERR] = LSECT5(PL1,PL2,TOL,CLIM) displays a
%         warning if the condition number of the matrix is greater than
%         a limit, CLIM.  This is usually due to parallel or nearly
%         parallel lines.  The default condition number limit is 1e+8.
%
%         NOTES:  1.  All the lines in the two piecewise linear lines
%                 are assumed not to be parallel to each other.
%
%                 2.  The piecewise linear lines must have point
%                 coordinates in two (2) columns (X,Y).
%
%                 3.  The M-files lsect3.m and lsect4.m must be in the
%                 current path or directory.
%
%                 4.  See lsect.m and lsect2.m for two-dimensional
%                 line (2-D) intersections.
%
%                 5.  See lisect.m and lisect2.m for three-dimensional
%                 line (3-D) intersections with infinite lines.  See
%                 lisect3.m for the intersection between two closed
%                 lines.
%
%         22-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error([' *** ERROR in  LSECT5:  LSECT5 requires two input', ...
         ' arguments.']);
end
%
if (nargin<3)||(isempty(tol))||(tol<=0)
  tol = 1e-8;
end
%
if (nargin<5)||(isempty(clim))||(clim<=1)
  clim = 1e+8;
end
%
% Check Points and Piecewise Linear Line (PLL)
%
[n1,l1] = size(pl1);
[n2,l2] = size(pl2);
%
if (n1<2)||(l1~=2)||(n2<2)||(l2~=2)
  error([' *** ERROR in LSECT5:  Error in input  piecewise linear ', ...
         'lines.']);
end
%
% Initialize Arrays
%
nl = n1-1;              % Number of lines in piecewise linear line (PL1)
%
ips = NaN(nl,2);
ts = NaN(nl,2);
idc = NaN(nl,2);
ierr = true;
%
% Loop through First Piecewise Linear Line
%
for k = 1:nl
%
   l = k+1;             % Next point
%
   p1 = pl1(k,:)';      % Point 1 on line segment
   p2 = pl1(l,:)';      % Point 2 on line segment
%
% Check for Intersection
%
   [ip,t,id2,ier] = lsect4(p1,p2,pl2,tol,clim);
   if ~isempty(ip)&&~ier
     ips(k,:) = ip';
     ts(k,:) = t';
     ierr = false;
     idc(k,:) = [k id2];
   end
%
end
idv = ~isnan(idc(:,1));
idc = idc(idv,:);
ips = ips(idv,:);
ts = ts(idv,:);
%
return