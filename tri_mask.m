function mask = tri_mask(pts,npx)
%TRI_MASK  Creates an image mask based on the two-dimensional (2-D)
%          coordinates of three points which form a triangle.  The
%          region of interest (ROI) is defined as the space within the
%          triangle
%          line.
%
%          MASK = TRI_MASK(PTS,NPX) given the two-dimensional
%          coordinates of three points defining a triangular region of
%          interest with X coordinates in the first column and Y in the
%          second column in PTS and the size of the image for the mask
%          in array, NPX,  creates a logical array mask for the region
%          of interest, MASK.  Note:  If NPX has only one element, the
%          image is assumed to be square (symmetrical) (NPX by NPX).
%
%          NOTES:  1.  M-file in_tri2d.m must be in the current
%                  directory or path.
%
%          04-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in TRI_MASK:  Two inputs are required!');
end
%
[nr,nc] = size(pts);
if nr~=3&&nc~=2
  error([' *** ERROR in TRI_MASK:  Input coordinate array must have'
         ' three rows and two columns!']);
end
%
ndim = size(npx(:),1);
if ndim>2||ndim<1||any(npx<1)
  error(' *** ERROR in TRI_MASK:  Incorrect image dimensions!');
end
%
if ndim==1
  npx = [npx; npx];
end
%
% Define Triangle Connectivity
%
tri = [1 2 3];
%
% Find Pixels within the Region of Interest
%
mask = false(npx(1)*npx(2),1);
%
minr = floor(min(pts));
maxr = ceil(max(pts));
idx = minr(:,1):maxr(:,1);
idy = minr(:,2):maxr(:,2);
[xg,yg] = meshgrid(idx,idy);
xym = [xg(:) yg(:)];
in_roi = in_tri2d(tri,pts,xym);
%
idr = sub2ind([npx(1) npx(2)],xym(:,2),xym(:,1));
idr = idr(in_roi);
%
mask(idr) = true;
%
return