function [tri,nt] = mk_tri4_2dc(dat,tol,iplt)
%MK_TRI4_2DC Makes a triangular mesh by using boundary line data from
%            a digitized MRI slice.  The first line is assumed to be
%            cartilage and the second line is assumed to be bone.
%
%        [TRI,NT] = MK_TRI4_2DC(DAT) given a cell array containing
%        two (2) columns matrices with boundary line coordinate point
%        data, DAT, returns the three (3) column triangle connectivity
%        matrix, TRI.  The number of returned triangles, NT, may also
%        be returned.
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
%                3.  The arclength along each line is used to determine
%                the triangulation.
%
%                4.   M-file lsect2.m must be in the current directory
%                or path.
%
%        16-Mar-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  iplt = false;
end
%
if nargin<2
  tol = 0.1;
end
if isempty(tol)
  tol = 0.1;
end
%
if (nargin<1)
  error(' *** ERROR in mk_tri4_2dc:  No input data!');
end
%
% Get Arc Lengths
%
dat = dat(:);
nslice = size(dat,1);
slen = cell(nslice,1);
npts = zeros(nslice,1);
rpt1 = zeros(nslice,2);
vec1 = zeros(nslice,2);
irev = false;
for k = 1:nslice
   xy = dat{k};
   vec = xy(end,:)-xy(1,:);
   vec1(k,:) = vec./norm(vec);
%
% Check for Slices with a Reverse Digitization
%
   if k>1
     dotp = vec1(k-1,:)*vec1(k,:)';
     if dotp<tol
       irev = true;
       xy = flipud(xy);
       vec = xy(end,:)-xy(1,:);
       vec1(k,:) = vec./norm(vec);
       dotp2 = vec1(k-1,:)*vec1(k,:)';
       if dotp2<dotp    % Revert back to original ordering
         warning([' *** WARNING in mk_tri4_2dc:  Ordering of points', ...
                  ' in the slices may not be in the same direction!']);
         irev = false;
         xy = flipud(xy);
         vec = xy(end,:)-xy(1,:);
         vec1(k,:) = vec./norm(vec);
       end
     else
       irev = false;
     end
   end
%
   rpt1(k,:) = xy(1,:);
   npts(k) = size(xy,1);
   dd = diff(xy);
   if k==1
     de = dd([1,npts(k)-1],:);         % Slopes at ends of top line
   end
   dlen = sqrt(sum(dd.*dd,2));
   slen{k} = [0; cumsum(dlen)];
   if irev
     slen{k} = flipud(slen{k});
   end
end
%
n = [0; cumsum(npts)];
tri = [];
slx = zeros(nslice-1,1);
for k = 2:nslice
%
% Slice Separations and Offsets
%
   ds = rpt1(k,:)-rpt1(k-1,:);
   offst = ds*vec1(k-1,:)';
   slx(k-1) = norm(ds)-norm(offst);
%
% Ends of Top (Cartilage) Line
%    mp = -de(:,1)./de(:,2);             % 90 degrees
   mp(2,1) = (de(2,1)+de(2,2))./(de(2,1)-de(2,2));    % 45 degrees
   mp(1,1) = (de(1,2)-de(1,1))./(de(1,1)+de(1,2));    % 45 degrees
   xp = dat{k-1}([1;npts(k-1)],1);
   yp = dat{k-1}([1;npts(k-1)],2);
   bp = yp-mp.*xp;
%
   xy = dat{k};
   if irev
     xy = flipud(xy);
   end
%
% Cut Off Extra Bone
% Fails if Bone is Not Longer Than Cartilage
%
   [~,~,id1] = lsect2(mp(1),bp(1),xy);
   if isempty(id1)||id1>npts(k)/4
     id1 = 1;
   end
   [~,~,id2] = lsect2(mp(2),bp(2),xy);
   if isempty(id2)||id2<npts(k)-npts(k)/4
     id2 = npts(k)-1;
   end
   idc = id1:id2+1;
   nptc = length(idc);
%
% Delaunay Triangulation
%
   xt = [zeros(npts(k-1),1); slx(k-1)*ones(nptc,1)];
   yt = [slen{k-1}-offst; slen{k}(idc)];
 %  xt = [zeros(npts(k-1),1); ones(npts(k),1)];
 %  yt = [(0:1/(npts(k-1)-1):1)'; (0:1/(npts(k)-1):1)'];
   tril = delaunay(xt,yt);
   nid = n(k-1)+1:n(k+1);
   tri = [tri; nid(tril)];
   if iplt
     h1 = figure;
     orient tall;
     ntril = size(tril,1);
     cla;
     plot(xt,yt,'k.');
     hold on;
     trimesh(tril,xt,yt);
     text(xt,yt,int2str((1:length(xt))'),'Color','b','FontSize',12);
     text(mean(xt(tril),2),mean(yt(tril),2),int2str((1:ntril)'), ...
          'Color','r','FontSize',12);
     h2 = figure;
     orient tall;
     nt = size(tri,1);
     cla;
     plot(xy(:,1),xy(:,2),'k.-');
     hold on;
     xyc = dat{k-1};
     plot(xyc(:,1),xyc(:,2),'r.-');
     xycp = [xyc(1,:); xyc(1,:)-de(1,:)];
     plot(xycp(:,1),xycp(:,2),'r:');
     xycp = [xyc(end,:); xyc(end,:)+de(2,:)];
     plot(xycp(:,1),xycp(:,2),'r:');
     xy1 = [xyc; dat{k}];
     xx = xy1(:,1);
     yy = xy1(:,2);
     xa = [min(xx)-2; max(xx)+2];
     yw = [min(yy)-2 max(yy)+2];
     ya = mp(1)*xa+bp(1);
     plot(xa,ya,'g-','Color',[0 0.7 0]);
     ya = mp(2)*xa+bp(2);
     plot(xa,ya,'g-','Color',[0 0.7 0]);
     xp = reshape(xy1(tri,1),nt,3)';
     yp = reshape(xy1(tri,2),nt,3)';
     xp = repmat(mean(xp),3,1)+0.75*(xp-repmat(mean(xp),3,1));
     yp = repmat(mean(yp),3,1)+0.75*(yp-repmat(mean(yp),3,1));
     patch(xp,yp,[0.0638 0.7446 0.7292]);
     text(xx,yy,int2str((1:length(xx))'),'Color','k', ...
          'FontSize',12);
     text(mean(xx(tril),2),mean(yy(tril),2),int2str((1:nt)'), ...
          'Color','b','FontSize',12);
     axlim = axis;
     axis equal;
     axis([axlim(1:2) yw]);
     pause;
     close(h1,h2);
   end
%
end
%
nt = size(tri,1);
%
return