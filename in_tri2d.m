function [in,intri] = in_tri2d(tri,xy,pts);
%IN_TRI2D Finds whether points are within the boundaries of 2-D
%         triangles.
%
%         IN = IN_TRI2D(TRI,XY,PTS) given a three (3) column triangle
%         connectivity matrix, TRI, a two (2) column X and Y coordinate
%         matrix of nodes within TRI, XY, and a two (2) column X and Y 
%         coordinate matrix of points to test whether they are within
%         the triangles, PTS, returns a logical vector with the same
%         number of rows as PTS, IN, that is true if the corresponding
%         point is within the triangles or false if not within the
%         triangles.
%
%         [IN,INTRI] = IN_TRI2D(TRI,XY,PTS) returns a column vector with
%         the index into TRI of the triangle that contains the
%         corresponding point.  If the point is not within any triangle
%         a zero is returned.
%
%         NOTES:  1.  The triangles are assumed to not overlap.
%
%         06-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  error([' *** ERROR in IN_TRI2D:  A triangle connectivity matrix' ...
         ' and two sets of XY coordinates are required as inputs!']);
end
%
[nr,nc] = size(tri);
if nc~=3
  error([' *** ERROR in IN_TRI2D:  The triangle connectivity matrix' ...
         ' must have three columns!']);
end
[nr2,nc2] = size(xy);
[npts,nc3] = size(pts);
if nc2~=2||nc3~=2
  error([' *** ERROR in IN_TRI2D:  The coordinate matrices' ...
         ' must have two columns!']);
end
%
% Initialize Variables for the Loop
%
x = xy(:,1);
y = xy(:,2);
%
in = false(npts,1);
intri = zeros(npts,1);
%
% Calculate the Divisor
%
del = (x(tri(:,2))-x(tri(:,1))).*(y(tri(:,3))-y(tri(:,1))) ...
      -(x(tri(:,3))-x(tri(:,1))).*(y(tri(:,2))-y(tri(:,1)));
%
% Loop through the Points
%
for k = 1:npts
   xi = pts(k,1);
   yi = pts(k,2);
%
% Get Point into Triangle Coordinates
%
   w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi)- ...
            (x(tri(:,1))-xi).*(y(tri(:,3))-yi))./del;
   w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi)- ...
            (x(tri(:,3))-xi).*(y(tri(:,2))-yi))./del;
   w12 = sum(w(:,1:2),2);
%
% Check to See if Within Triangles
%
   lidx = w(:,1)>=0&w(:,1)<=1&w(:,2)>=0&w(:,2)<=1&w12<=1;
   in(k) = any(lidx);
   if nargout>1&&in(k)
     intri(k) = find(lidx);
   end
end
%
return