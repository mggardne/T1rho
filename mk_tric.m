function [tri,xyc,nt] = mk_tric(xy);
%MK_TRIC Makes a triangular mesh using 2D points that form a closed
%        loop.
%
%        [TRI,XYC] = MK_TRIC(XY) given a matrix with two (2) columns
%        with coordinate point data, XY, returns the three (3) columns
%        triangle connectivity matrix, TRI and new coordinates that
%        includes the center of the loop, XYC.
%
%        [TRI,XYC,NT] = The number of returned triangles, NT, may also
%        be returned.
%
%        NOTES:  1.  X coordinates should be in column 1 and Y
%                coordinates should be in column 2.
%
%                2.  The coordinates should be ordered in the same
%                direction around the loop.
%
%                3.  All of the coordinates should be unique.  The
%                function will close the loop between the last and
%                first point.
%
%        16-Mar-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in MK_TRIC:  No input data!');
end
%
% Get Center of Loop and New Coordinates
%
xyc = [xy; mean(xy)];
%
% Create Triangles
%
nt = size(xy,1);        % Number of triangles equal the number of original points
nptsp1 = nt+1;
nptsm1 = nt-1;
%
tri = [repmat(nptsp1,nptsm1,1) [1:nptsm1]' [2:nt]'];
tri = [tri; [nptsp1 nt 1]];
%
return