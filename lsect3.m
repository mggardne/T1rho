function [ip,t,ierr] = lsect3(p1,v1,p2,v2,tol,clim)
%LSECT3   Finds the intersection of two closed two-dimensional (2-D)
%         lines.
%
%         IP = LSECT3(P1,V1,P2,V2) finds the intersection of two lines
%         defined by a point (P1) and a vector (V1) for the first line
%         and by a point (P2) and a vector (V2) for the second line.
%         The X and Y coordinates of the intersection are returned in
%         IP.
%
%         [IP,T] = LSECT3(P1,V1,P2,V2) returns the distance along the
%         lines to the intersection.  T(1) is the normalized (0 to 1)
%         distance along the first line to the intersection.  T(2) is
%         the distance along the second line to the intersection. 
%
%         [IP,T,IERR] = LSECT3(P1,V1,P2,V2,TOL) sets IERR to true if
%         no intersection is found within tolerance (TOL).  Default
%         tolerance is 1e-8.
%
%         [IP,T,IERR] = LSECT3(P1,V1,P2,V2,TOL,CLIM) sets IERR to true
%         if the condition number of the matrix is greater than a
%         limit, CLIM.  This is usually due to parallel or nearly
%         parallel lines.  The default condition number limit is 1e+8.
%
%         NOTES:  1.  The lines are assumed to be not parallel.
%
%                 2.  All input points and vectors must be of length
%                 two (2).
%
%                 3.  See lsect.m and lsect2.m for two-dimensional
%                 line (2-D) intersections with infinite lines.
%
%                 4.  See lisect.m and lisect2.m for three-dimensional
%                 line (3-D) intersections with infinite lines.
%
%         19-Nov-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<4
  error([' *** ERROR in LSECT3:  LSECT3 requires four input', ...
         ' arguments.']);
end
%
if (nargin<5)||(isempty(tol))||(tol<=0)
  tol = 1e-8;
end
%
if (nargin<6)||(isempty(clim))||(clim<=1)
  clim = 1e+8;
end
%
% Check Vectors
%
p1 = p1(:);
v1 = v1(:);
p2 = p2(:);
v2 = v2(:);
n1 = size(p1,1);
n2 = size(v1,1);
n3 = size(p2,1);
n4 = size(v2,1);
%
if (n1~=2)||(n2~=2)||(n3~=2)||(n4~=2)
  error([' *** ERROR in LSECT3:  All input points and vectors' ...
         ' must be of length two (2).']);
end
%
% Set Error Code and Set Outputs to NaNs
%
ierr = true;
ip = NaN(2,1);
t = NaN(2,1);
%
% Intersection Equations
%
A = [v1(1:2) -v2(1:2)];                % Solve for intersection
b = p2(1:2)-p1(1:2);
%
if rank(A)~=2
  disp(' *** WARNING in LSECT3:  Lines are parallel!');
  return
end
%
% Check Condition Number
%
cn = cond(A);
if cn>clim
  disp([' *** WARNING in LSECT3:  Lines are parallel' ...
            ' within the condition number limit.']);
  return
end
%
% Solve for Intersection
%
tv = A\b;
%
xy1 = tv(1)*v1+p1;                     % Coordinates of intersection
xy2 = tv(2)*v2+p2;
%
chk = norm(xy1-xy2);                   % Check intersection
if chk>tol
  disp([' *** WARNING in LSECT3:  Lines do not intersect' ...
            ' within tolerance.']);
  return
else
  if all(tv>=0)&&all(tv<=1)
    ip = xy1;
    t = tv;
    ierr = false;
  end
end
%
return