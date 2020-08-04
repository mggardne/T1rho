%#######################################################################
%
%                     * T1rho SLice PLoT Program *
%
%          M-File which reads knee MRI T1rho data from MAT files and
%     OsiriX digitized cartilage regions of interests (ROIs) from CSV
%     files to generate plots of the cartilage T1rho values and masks
%     for the cartilage and femoral condyle bone data.  Masks, ROIs and
%     slice information are saved into the original MAT file.
%
%     NOTES:  1.  The T1rho MAT files must be in the current directory.
%
%             2.  M-files cr_mask.m, in_tri2d.m, lsect2.m, meshbnd3.m,
%             mk_tri4_2dc.m, mk_tric.m and rd_roi5.m must be in the
%             current directory or path.
%
%     03-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Pick MRI Series to Analyze
%
mnams = dir('T1rho*.mat');
mnams = {mnams.name}';
%
idm = menu('Pick a MAT File to Analyze',mnams);
mnam = mnams{idm};
%
fs = extractAfter(mnam,'_');
fs = fs(1:end-4);
%
% Load T1rho Data
%
load(mnam,'npx','nslt','T1rhonls');
nsl = size(T1rhonls,3); % Number of slices
%
% Get ROIs
%
rdir = fullfile('JULYT1Rho',fs);
rnams = dir(fullfile(rdir,'*.csv'));
rnams = {rnams.name}';
nrfiles = size(rnams,1);
%
% Loop through ROI Files
%
rois = struct;
%
for k = 1:nrfiles
%
   rnam = fullfile(rdir,rnams{k});
%
   rois(k).name = rnams{k};
   rois(k).roi = rd_roi5(rnam,true);
   rois(k).slice = [rois(k).roi.imageno]';
%
end
%
% Get Slices with ROIs
%
rsl = {rois.slice}';
rsl = cell2mat(rsl);
rsl = unique(rsl);
idf = (1:nslt:nslt*nsl)';              % Spin lock time = 0 ms slices
%
% Check for Slices Greater Than Spin Lock Time = 0 ms and Set to
% Slice Numbers with Associated Spin Lock Time = 0 ms
%
[ide,idxe] = setdiff(rsl,idf);
if ~isempty(ide)
  ne = size(ide,1);
  for k = 1:ne
%
     idfl = find(idf<ide(k));
     idfl = idfl(end);
     idfl = idf(idfl);
     rsl(idxe(k)) = idfl;
%
% Correct ROIS to Point to Correct Slices
%
     for l = 1:nrfiles
        le = rois(l).slice==ide(k);
        if any(le)
          rois(l).slice(le) = idfl;
          rois(l).roi(le).imageno = idfl;
        end
     end
%
  end
%
  rsl = unique(rsl);    % Get unique number of slices
%
end
%
% Convert From Spin Lock Time = 0 ms to T1rho Slices
%
nrsl = size(rsl,1);
[~,~,irsl] = intersect(rsl,idf);       % T1rho slices with ROIs
%
% Get Index to Bone and Cartilage Regions 
%
legds = strrep(rnams,'.csv','');
rnamsc = char(legds);
bones = rnamsc(:,1);    % Initial letters of femur, patella and tibia
bc = rnamsc(:,4);       % Initial letters of bone or cartilage
%
% Setup Color Maps
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = zeros(128,3);    % Jet color map for cartilage
jmap(13:90,:) = jet(78);% Jet colors in scaled range 128*[0 10 70 100]./100 
jmap(1:12,:) = repmat(jmap(13,:),12,1);          % Set lower values to blue
jmap(91:128,:) = repmat(jmap(90,:),128-90,1);    % Set higher values to red
cmap = [gmap; jmap];
%
bmap = gray(256);       % Color map for bone area
bmap(1,:) = [0.8 0 0];              % Bone is red
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
lt = ['b.-'; 'g.-'; 'r.-'; 'c.-'; 'm.-'; 'y.-' ];  % Line color and type
pnam1 = [fs '_T1rho_ROIs1.ps'];        % ROI lines print file name
pnam2 = [fs '_T1rho_ROIs2.ps'];        % ROI areas print file name
hnam = [fs '_Histograms.ps'];          % Histograms print file name
bnam = [fs '_BoneROIs.ps'];            % Bone ROIs print file name
%
f = cell(2,nrsl);       % Femur coordinates (1 - cartilage, 2 - bone)
p = cell(2,nrsl);       % Patella coordinates (1 - cartilage, 2 - bone)
t = cell(2,nrsl);       % Tibia coordinates (1 - cartilage, 2 - bone)
%
ibone = false(nrsl,3);  % 1 - femur, 2- patella and 3 - tibia
maskf = false(npx*npx,nrsl);           % Mask for femoral cartilage
maskp = false(npx*npx,nrsl);           % Mask for patellar cartilage
maskt = false(npx*npx,nrsl);           % Mask for tibial cartilage
%
maskb = false(npx*npx,nrsl);           % Mask for femoral condyle bone
%
for k = 1:nrsl
%
% Plot T1rho Slice Images
%
   slk = rsl(k);        % Slice number in spin lock images
   sl = irsl(k);        % Slice number in T1rho images
%
   figure;
   orient landscape;
   imagesc(T1rhonls(:,:,sl),[0 100]);
   colormap gray;
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
   hold on;
%
   lh = gobjects(nrfiles,1);           % Line graphic handles
   idl = false(nrfiles,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
   for l = 1:nrfiles
      idxr = rois(l).slice==slk;
      if any(idxr)
        dat = cell2mat(rois(l).roi(idxr).data);
        b = bones(l);
        m = double(strcmpi(bc(l),'B'))+1;
        if strcmpi(b,'F')
          f{m,k} = dat; % Femur
          ibone(k,1) = true;
        elseif strcmpi(b,'P')
          p{m,k} = dat; % Patella
          ibone(k,2) = true;
        else
          t{m,k} = dat; % Tibia
          ibone(k,3) = true;
        end
        lh(l) = plot(dat(:,1),dat(:,2),lt(l,:));
        idl(l) = true;
      end
   end
%
% Add Legends and Print Slice Plots
%
   legend(lh(idl),legds(idl));
   if k==1
     print('-dpsc2','-r600','-fillpage',pnam1);
   else
     print('-dpsc2','-r600','-fillpage','-append',pnam1);
   end
%
% Create Logical Masks for the Cartilage on this Slice
%
   if ibone(k,1)
     maskf(:,k) = cr_mask(f(:,k),npx);
   end
   if ibone(k,2)
     maskp(:,k) = cr_mask(p(:,k),npx);
   end
   if ibone(k,3)
     maskt(:,k) = cr_mask(t(:,k),npx);
   end
%
% Plot ROIs
%
   mask = maskf(:,k)|maskp(:,k)|maskt(:,k);      % All cartilage mask
   img = T1rhonls(:,:,sl);
   idmx = img>100;
   img(idmx) = 100;
   idmn = img<0;
   img(idmn) = 0;
   imgb = img;          % Image data for bone ROI (below)
   img(~mask) = img(~mask)-100.1;
%
   figure;
   orient landscape;
   imagesc(img,[-100 100]);
   colormap(cmap);
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
   hb = colorbar;
   set(hb,'Limits',[10 70]);
%
   if k==1
     print('-dpsc2','-r600','-fillpage',pnam2);
   else
     print('-dpsc2','-r600','-fillpage','-append',pnam2);
   end
%
% Plot Histograms of the Cartilage T1rho Values
%
   figure;
   orient landscape;
%
   histogram(img(mask),25,'FaceAlpha',1,'FaceColor',[0 0 0.8]);
   xlabel('T1rho Values (ms)','FontSize',12,'FontWeight','bold');
   ylabel('Frequency','FontSize',12,'FontWeight','bold');
   title([fs ' Slice ' int2str(sl) ' Histogram'],'FontSize',16, ...
         'FontWeight','bold');
%
   if k==1
     print('-dpsc2','-r600','-fillpage',hnam);
   else
     print('-dpsc2','-r600','-fillpage','-append',hnam);
   end
%
% Create Logical Masks for the Femoral Condyle Bone on this Slice
%
   datb = f{2,k};       % Femoral bone line
%
   [trib,xyb] = mk_tric(datb);
   o = meshbnd3(trib);  % Get boundary nodes
%
   minr = round(min(xyb));
   maxr = round(max(xyb));
   idx = minr(:,1):maxr(:,1);
   idy = minr(:,2):maxr(:,2);
   [xg,yg] = meshgrid(idx,idy);
   xym = [xg(:) yg(:)];
   in_b = in_tri2d(trib,xyb,xym);
%
   idb = sub2ind([npx npx],xym(:,2),xym(:,1));
   idb = idb(in_b);
%
   maskb(idb,k) = true;
%
% Plot Femoral Condyle Bone Outline and Area
%
   figure;
   orient landscape;
   imagesc(imgb,[0 100]);
   colormap gray;
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
   hold on;
   plot(xyb(o,1),xyb(o,2),'r.-');
%
   if k==1
     print('-dpsc2','-r600','-fillpage',bnam);
   else
     print('-dpsc2','-r600','-fillpage','-append',bnam);
   end
%
   figure;
   orient landscape;
   imgm = imgb;
   imgm(maskb(:,k)) = -0.5;
   imagesc(imgm,[-0.5 100]);
   colormap(bmap);
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
%
   print('-dpsc2','-r600','-fillpage','-append',bnam);
%
end
%
% Save Masks, ROIS and Slice Information into MAT File
%
save(mnam,'bc','bones','f','ibone','irsl','maskb','maskf','maskp', ...
     'maskt','nrfiles','p','rois','rsl','t','-append');
%
return