function [mask1,mask2] = cr_mask2(roic,npx,dist,scal,tol,iplt)
%CR_MASK2  Creates two image masks based on the two-dimensional (2-D)
%          coordinates of two lines.  The region of interests (ROIs)
%          are defined as the space between the first line and the
%          midline between the first and second lines and the space
%          between the midline and the second line.  The second line is
%          assumed to be slightly longer than the first line.
%
%          For creating masks for knee joint cartilage.  The first line
%          is usually cartilage and the second longer line is the
%          underlying bone.
%
%          [MASK1,MASK2] = CR_MASK2(ROIC,NPX) given the two-dimensional
%          coordinates of two lines defining a region of interest in a
%          two element cell array (one element for each line), ROIC,
%          and the size of the image for the mask in array, NPX,
%          creates two logical array mask for the region of interest,
%          MASK1 and MASK2.  Note:  If NPX has only one element, the
%          image is assumed to be square (symmetrical) (NPX by NPX).
%
%          [MASK1,MASK2] = CR_MASK2(ROIC,NPX,DIST,SCAL) midline points
%          must be within a distance DIST of the second line.  The
%          default value is Inf.  A one or two element scale SCAL is
%          used to scale the X and Y coordinates before comparison with
%          the distance DIST for midline points.  See midline.m.
%
%          NOTES:  1.  M-files in_tri2d.m, lsect2.m, lsect2a.m,
%                  midline.m, mk2_tri_2d.m must be in the current
%                  directory or path.
%
%          21-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in CR_MASK2:  Two inputs are required!');
end
%
roic = roic(:);
[nr,nc] = size(roic);
if nr~=2&&nc~=1
  error([' *** ERROR in CR_MASK2:  Input cell array must have'
         ' two elements!']);
end
%
ndim = size(npx(:),1);
if ndim>2||ndim<1
  error([' *** ERROR in CR_MASK2:  Incorrect number of image ', ...
         'dimensions!']);
end
%
if ndim==1
  npx = [npx; npx];
end
%
if nargin<3
  dist = Inf;            % No distance checking
end
%
if isempty(dist)
  dist = Inf;
end
%
if nargin<4
  scal = [1 1];         % No scaling
end
if isempty(scal)
  scal = [1 1];
end
scal = scal(:);
nr = size(scal,1);
if nr==1
  scal = [scal scal];   % Same scaling for X and Y
else
  scal = scal(1:2)';
end
%
if nargin<5
  tol = 0.1;
end
if isempty(tol)
  tol = 0.1;
end
%
if nargin<6
  iplt = false;
end
%
% Create Triangle Meshes for the Region of Interest
%
[tri1,tri2,xy] = mk2_tri_2d(roic,dist,scal,tol,iplt);
%
% Find Pixels within the Region of Interest
%
mask1 = false(npx(1)*npx(2),1);
mask2 = false(npx(1)*npx(2),1);
%
minr = floor(min(xy));
if any(minr(:)==0)      % Trap for zeros
  minr(minr(:)==0) = 1;
end
maxr = ceil(max(xy));
idx = minr(:,1):maxr(:,1);
idy = minr(:,2):maxr(:,2);
[xg,yg] = meshgrid(idx,idy);
xym = [xg(:) yg(:)];
%
idr = sub2ind([npx(1) npx(2)],xym(:,2),xym(:,1));
%
in_roi = in_tri2d(tri1,xy,xym);
idr1 = idr(in_roi);
mask1(idr1) = true;
%
in_roi = in_tri2d(tri2,xy,xym);
idr2 = idr(in_roi);
mask2(idr2) = true;
%
return