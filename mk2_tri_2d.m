function [tri1,tri2,xy,xym] = mk2_tri_2d(dat,dist,scal,tol,iplt)
%MK2_TRI_2D Makes two triangular meshes by using boundary line data from
%           a two-dimensional digitized MRI slice.  The first line is
%           assumed to be cartilage and the second line is assumed to
%           be bone.  The region is divided into two regions using a
%           midline between the two lines.
%
%        [TRI1,TRI2,XY] = MK2_TRI_2D(DAT) given a cell array
%        containing two (2) columns matrices with boundary line
%        coordinate point data, DAT, returns two three (3) column
%        triangle connectivity matrices, TRI1 and TRI2 and X and Y
%        coordinates in two columns matrices XY1 and XY2.
%
%        [TRI1,TRI2,XY] = MK2_TRI_2D(DAT,DIST) midline points must
%        be within a distance DIST of the second line.  The default
%        value is Inf.  See midline.m.
%
%        [TRI1,TRI2,XY] = MK2_TRI_2D(DAT,DIST,SCAL) a one or two
%        element scale SCAL is used to scale the X and Y coordinates
%        before comparison with the distance DIST for midline points.
%        See midline.m.
%
%        NOTES:  1.  Each boundary coordinate data matrix must
%                correspond to one index into the cell array DAT.
%
%                2.  The coordinates should be ordered in the same
%                direction in each line.  The dot product of the
%                directions of the adjacent lines are used to check the
%                ordering direction and the ordering direction is
%                reversed if the dot product is negative.
%
%                3.  The first triangulation is between the first line
%                and the midline.  The second triangulation is between
%                the midline and the second line.
%
%                3.   M-files midline.m, lsect2.m and lsect2a.m must be
%                in the current directory or path.
%
%        20-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in mk2_tri_2d:  No input data!');
end
%
if nargin<2
  dist = Inf;            % No distance checking
end
%
if isempty(dist)
  dist = Inf;
end
%
if nargin<3
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
if nargin<4
  tol = 0.1;
end
if isempty(tol)
  tol = 0.1;
end
%
if nargin<5
  iplt = false;
end
%
dat = dat(:);
nslice = size(dat,1);
%
if nslice~=2
  error([' *** ERROR in mk2_tri_2d:  Input cell array must ', ...
         'have two elements containing 2D coordinates for two lines!']);
end
%
% Get First (Cartilage) Line and Slopes at the Ends of the Line
%
xy1 = dat{1};
npts1 = size(xy1,1);
%
vec1 = xy1(npts1,:)-xy1(1,:);          % Direction of line
vec1 = vec1./norm(vec1);
%
dd = diff(xy1);
de = dd([1,npts1-1],:);         % Slopes at ends of top line
%
% Ends of First (Cartilage) Line
%
% mp = -de(:,1)./de(:,2); % 90 degrees
mp(2,1) = (de(2,1)+de(2,2))./(de(2,1)-de(2,2));  % 45 degrees
mp(1,1) = (de(1,2)-de(1,1))./(de(1,1)+de(1,2));  % 45 degrees
xp = xy1([1; npts1],1);
yp = xy1([1; npts1],2);
bp = yp-mp.*xp;
%
% Get Second (Bone) Line
%
xy2 = dat{2};
npts2 = size(xy2,1);
%
vec2 = xy2(npts2,:)-xy2(1,:);
vec2 = vec2./norm(vec2);
%
% Check for Slices with a Reverse Digitization
%
dotp = vec1*vec2';
%
if dotp<tol
  xy2 = flipud(xy2);
  vec = xy2(npts2,:)-xy2(1,:);
  vec = vec./norm(vec);
  dotp2 = vec1*vec';
  if dotp2<dotp         % Revert back to original ordering
    warning([' *** WARNING in mk2_tri4_2d:  Ordering of points', ...
             ' in the slices may not be in the same direction!']);
    xy2 = flipud(xy2);
  end
end
%
% Cut Off Extra Points on Second (Bone) Line
% Fails if Second (Bone) is Not Longer Than First (Cartilage) Line
%
[~,~,id1] = lsect2(mp(1),bp(1),xy2);
if isempty(id1)||id1>npts2/4
  id1 = 1;
end
[~,~,id2] = lsect2(mp(2),bp(2),xy2);
if isempty(id2)||id2<npts2-npts2/4
  id2 = npts2-1;
end
idc = id1:id2+1;
nptc = length(idc);
xy2 = xy2(idc,:);
%
% Get Midline
%
xym = midline(xy1,xy2,dist,scal,iplt);
nptsm = size(xym,1);
%
xy = [xy1; xym; xy2];
%
% Delaunay Triangulation
%
n = [0; cumsum([npts1; nptsm; nptc])];
%
c1 = [(n(1)+1:n(2)-1)' (n(1)+2:n(2))'; n(2) n(3); (n(3):-1:n(2)+2)' ...
      (n(3)-1:-1:n(2)+1)'; n(2)+1 n(1)+1];       % Constraints
%
dt1 = delaunayTriangulation(xy,c1);
idin = isInterior(dt1);
tri1 = dt1(idin,:);
%
c2 = [(n(2)+1:n(3)-1)' (n(2)+2:n(3))'; n(3) n(4); (n(4):-1:n(3)+2)' ...
      (n(4)-1:-1:n(3)+1)'; n(3)+1 n(2)+1];       % Constraints
%
dt2 = delaunayTriangulation(xy,c2);
idin = isInterior(dt2);
tri2 = dt2(idin,:);
%
% Plot Triangulations?
%
if iplt
%
  h1 = figure;
  orient tall;
%
  nt1 = size(tri1,1);
  nt2 = size(tri2,1);
  xt = xy(:,1);
  yt = xy(:,2);
%
  plot(xt,yt,'k.','LineWidth',1,'MarkerSize',7);
  hold on;
  npts = size(xt,1);
  text(xt,yt,int2str((1:npts)'),'Color','k','FontSize',10);
%
  trimesh(tri1,xt,yt);
  text(mean(xt(tri1),2),mean(yt(tri1),2),int2str((1:nt1)'), ...
       'Color','r','FontSize',10);
%
  trimesh(tri2,xt,yt);
  text(mean(xt(tri2),2),mean(yt(tri2),2),int2str((1:nt2)'), ...
       'Color','b','FontSize',10);
%
  h2 = figure;
  orient tall;
  plot(xy2(:,1),xy2(:,2),'k.-','LineWidth',1,'MarkerSize',7);
  hold on;
  plot(xy1(:,1),xy1(:,2),'b.-','LineWidth',1,'MarkerSize',7);
  plot(xym(:,1),xym(:,2),'r.:','LineWidth',1,'MarkerSize',7);
  text(xt,yt,int2str((1:npts)'),'Color','k', ...
       'FontSize',10);
%
  xp = reshape(xy(tri1,1),nt1,3)';
  yp = reshape(xy(tri1,2),nt1,3)';
  xp = repmat(mean(xp),3,1)+0.75*(xp-repmat(mean(xp),3,1));
  yp = repmat(mean(yp),3,1)+0.75*(yp-repmat(mean(yp),3,1));
  patch(xp,yp,[1 0.7 0.7]);
  text(mean(xt(tri1),2),mean(yt(tri1),2),int2str((1:nt1)'), ...
       'Color','r','FontSize',10);
%
  xp = reshape(xy(tri2,1),nt2,3)';
  yp = reshape(xy(tri2,2),nt2,3)';
  xp = repmat(mean(xp),3,1)+0.75*(xp-repmat(mean(xp),3,1));
  yp = repmat(mean(yp),3,1)+0.75*(yp-repmat(mean(yp),3,1));
  patch(xp,yp,[0.3 0.7 1]);
  text(mean(xt(tri2),2),mean(yt(tri2),2),int2str((1:nt2)'), ...
       'Color','b','FontSize',10);
%
  axis equal;
  pause;
  close(h1,h2);
%
end
%
return