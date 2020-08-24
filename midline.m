function xymid = midline(xyc,xyb,tol,pspc,iplt)
%MIDLINE   Creates a midline between the two-dimensional (2-D)
%          coordinates of two lines.  The midline is based on
%          perpendiculars to the second line.  The second line is
%          assumed to be slightly longer than the first line.   The two
%          lines must be going in the same direction so that the two
%          ends of the lines line up.
%
%          For creating masks for knee joint cartilage.  The first line
%          is usually cartilage and the second longer line is the
%          underlying bone.
%
%          XYMID = MIDLINE(XYC,XYB) given the two-dimensional
%          coordinates of two lines defining a region of interest in
%          two column matrices XYC and XYB, the two-dimensional
%          coordinates of the midline are returned in the two column
%          matrix XYMID.
%
%          XYMID = MIDLINE(XYC,XYB,TOL) only finds midline points
%          within a tolerance of TOL of the points in the second line
%          XYB.
%
%          XYMID = MIDLINE(XYC,XYB,TOL,PSPC) scales the X and Y
%          coordinates before checking whether midline points are
%          within a tolerance of TOL of the points in the second line
%          XYB.  The first element of PSPC is used to scale X
%          and the second element is used to scale the Y coordinates.
%          Note:  If PSPC has only one element, the scaling in assumed
%          to be the same in X and Y.
%
%          NOTES:  1.  The two lines must be going in the same
%                  direction so that the two ends line up.
%
%                  2.  The midline is based on perpendiculars to the
%                  second line.
%
%                  3.  M-files lsect2a.m, lsect3.m, lsect4.m and
%                  lsect3.m, must be in the current directory or path.
%
%          20-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<2
  error(' *** ERROR in MIDLINE:  Two inputs are required!');
end
%
if nargin<3
  tol = Inf;            % No tolerance checking
end
%
if isempty(tol)
  tol = Inf;
end
%
tol = 2*tol;
tol = tol(1)*tol(1);    % Get square of twice tolerance for comparison
%
if nargin<4
  pspc = [1 1];         % No pixel scaling
end
%
if isempty(pspc)
  pspc = [1 1];
end
%
pspc = pspc(:);
nr = size(pspc,1);
%
if nr==1
  pspc = [pspc pspc];
else
  pspc = pspc(1:2)';
end
%
if nargin<5             % No plotting
  iplt = false;
end
%
% Get Slopes of the Second (Bone) Line
%
nb = size(xyb,1);       % Number of bone points
ni = nb-2;              % Number of interior points
%
dbb = diff(xyb);
dfb = dbb(1:nb-2,:);    % Forward slope
dbb = dbb(2:nb-1,:);    % Backward slope
dcb = (dfb+dbb)/2;      % Central slope
%
% Check Slopes
%
if any(abs(dcb(:,2))<1e-11)
  idx = abs(dcb(:,2))<1e-11;
  if any(abs(dbb(idx,2))>1e-11)
    idx0 = abs(dbb(idx,2))>1e-11;
    idx0 = idx(idx0);
    dcb(idx0,:) = dbb(idx0,2);         % Use backward slope
  end
  if any(abs(dfb(idx,2))>1e-11)
    idx0 = abs(dfb(idx,2))>1e-11;
    idx0 = idx(idx0);
    dcb(idx0,:) = dfb(idx0,2);         % Use forward slope
  end
  if any(abs(dcb(:,2))<1e-11)
    idx = abs(dcb(:,2))<1e-11;
    dcb(idx,:) = [1e+4 1];             % Use a steep slope
  end
end
%
% Get Perpendiculars
%
idi = 2:nb-1;           % Index to interior points
mb = -dcb(:,1)./dcb(:,2);    % Get perpendicular slopes
xb = xyb(idi,1);
yb = xyb(idi,2);
bb = yb-mb.*xb;         % Y-intercept for lines
%
% Plot Points
%
if iplt
  hf = figure;
  plot(xyc(:,1),xyc(:,2),'b.-');
  axis equal;
  hold on;
  plot(xyb(:,1),xyb(:,2),'k.-');
  %
  xp = [xb-15 xb+15]';
  yp = repmat(mb,1,2)'.*xp+repmat(bb,1,2)';
  plot(xp,yp,'g-');
  axis([min(xyc(:,1))-12 max(xyc(:,1))+12 min(xyc(:,2))-12 ...
        max(xyc(:,2))+12]);
end
%
% Find Intersections within Tolerance
%
xyi = zeros(nb-2,2);
for k = 1:ni
   xy = lsect2a(mb(k),bb(k),xyc);
   nxy = size(xy,1);
   d = xy-repmat(xyb(k,:),nxy,1);
   d = d.*repmat(pspc,nxy,1);
   d = sum(d.*d,2);
   [dmin,idmin] = min(d);              % Intersection closest to bone point
   if dmin<tol          % Within tolerance of bone point
     xyi(k,:) = xy(idmin,:);
   else
     xyi(k,:) = [NaN NaN];
   end
end
%
xymid = (xyi+xyb(idi,:))./2;           % Get midpoint
%
% Check for NaNs
%
bgpt = 1;
id1 = 0;
endpt = nb;
id2 = ni+1;
%
n1 = floor(ni/2);
%
idn = isnan(xymid(:,1));
idn1 = idn(1:n1);
idn2 = idn(n1+1:ni);
%
if any(idn1)
  id1 = find(idn1,1,'last');
  bgpt = idi(id1);
end
%
if any(idn2)
  id2 = find(idn2,1)+n1;
  endpt = idi(id2);
end
%
% Get End Points
%
xymid = [(xyb(bgpt,:)+xyc(1,:))/2; xymid(id1+1:id2-1,:); ...
         (xyb(endpt,:)+xyc(end,:))/2];
%
% Check that Midline Does Not Intersect First or Second Lines
%
[~,~,idi] = lsect5(xymid,xyc);         % Check against first line
if ~isempty(idi)
  nm = size(xymid,1);
  n1 = round(nm/2);     % Intersections should be near ends
  id1 = find(idi<n1,1,'last');
  if ~isempty(id1)
    id1 = idi(id1)+1;
    xymid = xymid(id1:nm,:);
  end
  id2 = find(idi>n1,1);
  if ~isempty(id2)
    id2 = idi(id2);
    xymid = xymid(1:id2,:);
  end
end
%
[~,~,idi] = lsect5(xymid,xyb);         % Check against second line
if ~isempty(idi)
  nm = size(xymid,1);
  n1 = round(nm/2);     % Intersections should be near ends
  id1 = find(idi<n1,1,'last');
  if ~isempty(id1)
    id1 = idi(id1)+1;
    xymid = xymid(id1:nm,:);
  end
  id2 = find(idi>n1,1);
  if ~isempty(id2)
    id2 = idi(id2);
    xymid = xymid(1:id2,:);
  end
end
%
% Finish Plot
%
if iplt
  plot(xymid(:,1),xymid(:,2),'ro:','MarkerSize',3, ...
       'MarkerFaceColor','r');
  pause;
  close(hf);
end
%
return