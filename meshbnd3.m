function out = meshbnd3(tri);
%MESHBND3 Finds and orders a list of boundary nodes so the nodes form a
%         closed connected boundary around the outside edges of a
%         triangular mesh.  See MESHBND2.M for a similar method using
%         the boundary nodes if they are known.
%
%         OUT = MESHBND3(TRI) given a three (3) column triangle
%         connectivity matrix, TRI, returns an ordered list of boundary
%         nodes, OUT, that forms a closed connected boundary of the
%         outside edges of the triangular mesh.
%
%         NOTES:  1.  The node list starts and ends with the lowest
%                 node ID.
%
%                 2.  The ordering can be either clockwise or counter-
%                 clockwise.
%
%                 3.  The node list will connect any holes in the mesh
%                 to other holes and the outside boundary.  Best for
%                 continuous surfaces without holes.
%
%                 4.  See MESHBND2.M for a similar method using the
%                 boundary nodes if they are known.  This is useful for
%                 limiting the outside to just the input boundary nodes
%                 (e.g. around a hole in the surface or just the
%                 exterior edges of a surface that includes a hole).
%
%                 5.  The M-file nod2tri.m must be in the current path
%                 or directory.
%
%         30-May-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<1
  error([' *** ERROR in MESHBND3:  A triangle connectivity matrix' ...
         ' is required as an input!']);
end
%
[nr,nc] = size(tri);
if nc~=3
  error([' *** ERROR in MESHBND3:  Triangle connectivity matrix' ...
         ' must have three columns!']);
end
%
% Get Sides of the Triangles
%
sides = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
sides = sort(sides,2);
ns = size(sides,1);
%
% Find Duplicate Edges (Duplicate Edge == Interior Edge)
%
sides = sortrows(sides);
ds = diff(sides);
ds = sum(abs(ds),2);
idup = find(ds==0);
%
% Delete Duplicate Edges
%
if ~isempty(idup)
  idup = [idup; idup+1];
  idx = true(ns,1);
  idx(idup) = false;
  sides = sides(idx,:);
  ns = size(sides,1);
end
%
% Get Lowest Node ID at the Start and End of the Outline
%
out = sides([1 3:ns 2]',:);
out(ns,:) = fliplr(out(ns,:));
%
% Order the Node IDs so that the node ID in column 2 of the row above
% the current row matches the node ID in column 1 of the current row
%
nr = ns-1;              % Final row is already set
for k = 2:nr;
   nn = out(k-1,2);
   idx = k:nr;          % Only look below the current row
   [i,j] = find(out(idx,:)==nn);       % Find match node ID
   if j==2
     out([k; idx(i)],:) = fliplr(out([idx(i); k],:));
   else
     out([k; idx(i)],:) = out([idx(i); k],:);
   end
end
%
out = [out(:,1); out(1)];
%
return