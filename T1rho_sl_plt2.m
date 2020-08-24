%#######################################################################
%
%                    * T1rho SLice PLoT 2 Program *
%
%          M-File which reads knee MRI T1rho data from MAT files and
%     OsiriX digitized cartilage regions of interests (ROIs) from CSV
%     files to generate plots of the cartilage T1rho values and masks
%     for the cartilage and femoral condyle bone data.  Masks, ROIs and
%     slice information are saved into the original MAT file.  N, mean,
%     standard deviation (SD) and signal-to-noise ratio (mean/SD) are
%     output for individual slices by bone and by layer (superficial or
%     deep) for the cartilage T1rho values within a valid range
%     (0<T1rho<100) are output to a MS-Excel spreadsheet.
%
%     NOTES:  1.  The T1rho MAT files must be in the current directory.
%             MAT file dicom_lst.mat must be in the current directory.
%             See dicom_lst.m.
%
%             2.  M-files cr_mask2.m, in_tri2d.m, lsect2.m, lsect2a.m,
%             lsect3.m, lsect4.m, lsect5.m, meshbnd3.m, midline.m,
%             mk2_tri_2d.m, mk_tric.m and rd_roi5.m must be in the
%             current directory or path.
%
%             3.  ROI CSV files must be in subdirectories under the
%             directory "\August T1 Segmentations".
%
%     19-Aug-2020 * Mack Gardner-Morse
%
%     21-Aug-2020 * Mack Gardner-Morse * Divided ROIs in the middle to
%                                        make two layers.
%
%     24-Aug-2020 * Mack Gardner-Morse * Added output of T1rho measures
%                                        to MS-Excel spreadsheet.
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
% Load T1rho Data
%
load(mnam,'fs','npx','nslt','sn','sn1','T1rhonls');
nsl = size(T1rhonls,3); % Number of slices
%
% Get Pixel Spacing
%
load('dicom_lst.mat','pspc');
idx = sn==sn1;
scal = pspc(idx,:);
dist = 7.5;             % Maximum distance to midline
%
% Get ROIs
%
rdir = 'August T1 Segmentations';
rsdirs = dir(rdir);
idd = [rsdirs.isdir]';
rsdirs = {rsdirs(idd).name}';
%
idd = startsWith(rsdirs,'.');          % Current and parent directories
rsdirs = rsdirs(~idd);  % Don't use current and parent directories
%
idr = menu('Pick Corresponding ROI Series Directory to Analyze',rsdirs);
rsdir = rsdirs{idr};
rdir = fullfile(rdir,rsdir);
%
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
idb = contains(rnams,'_')+4;           % Underscore between name parts?
bc = rnamsc(:,idb);     % Initial letters of bone or cartilage
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
roitxt = ['\fontsize{15}Superficial Layer is \color{blue}Blue', ...
          '\color{black} and Deep Layer is \color{red}Red'];
pnam3 = [fs '_T1rho_ROIs3.ps'];        % ROI values print file name
hnam = [fs '_Histograms.ps'];          % Histograms print file name
bnam = [fs '_BoneROIs.ps'];            % Bone ROIs print file name
%
f = cell(2,nrsl);       % Femur coordinates (1 - cartilage, 2 - bone)
p = cell(2,nrsl);       % Patella coordinates (1 - cartilage, 2 - bone)
t = cell(2,nrsl);       % Tibia coordinates (1 - cartilage, 2 - bone)
%
ibone = false(nrsl,3);  % 1 - femur, 2- patella and 3 - tibia
maskf = false(npx*npx,2,nrsl);         % Mask for femoral cartilage
maskp = false(npx*npx,2,nrsl);         % Mask for patellar cartilage
maskt = false(npx*npx,2,nrsl);         % Mask for tibial cartilage
%
maskb = false(npx*npx,nrsl);           % Mask for femoral condyle bone
%
% MS-Excel File
%
xnam = ['T1rho_sl_' fs '.xlsx'];
datas = cell(2,3,nrsl); % Data for all slices:  2 = layers and 3 = bones
slnums = zeros(2,4,nrsl);
bone_idx = zeros(2,4,nrsl);
nums = zeros(2,4,nrsl);
T1rho_mn = zeros(2,4,nrsl,'single');
T1rho_sd = zeros(2,4,nrsl,'single');
T1rho_snr = zeros(2,4,nrsl,'single');
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
   legend(lh(idl),legds(idl),'Interpreter','none');
   if k==1
     print('-dpsc2','-r600','-fillpage',pnam1);
   else
     print('-dpsc2','-r600','-fillpage','-append',pnam1);
   end
%
% Create Logical Masks for the Cartilage on this Slice
%
   if ibone(k,1)
     [maskf(:,1,k),maskf(:,2,k)] = cr_mask2(f(:,k),npx,dist,scal);
   end
   if ibone(k,2)
     [maskp(:,1,k),maskp(:,2,k)] = cr_mask2(p(:,k),npx,dist,scal);
   end
   if ibone(k,3)
     [maskt(:,1,k),maskt(:,2,k)] = cr_mask2(t(:,k),npx,dist,scal);
   end
%
% Plot ROIs
%
   mask1 = maskf(:,1,k)|maskp(:,1,k)|maskt(:,1,k);    % All cartilage mask
   mask2 = maskf(:,2,k)|maskp(:,2,k)|maskt(:,2,k);    % All cartilage mask
   img1 = T1rhonls(:,:,sl);
   idmx = img1>100;
   img1(idmx) = 100;
   idmn = img1<0;
   img1(idmn) = 0;
   imgb = img1;         % Image data for bone ROI (below)
   dats = img1;
   idmsk = ~mask1&~mask2;
   img1(idmsk) = img1(idmsk)-100.1;
   img2 = img1;
   img1(mask1) = 20;
   img1(mask2) = 60;
%
   figure;
   orient landscape;
   imagesc(img1,[-100 100]);
   colormap(cmap);
   axis image;
   axis off;
   title({[fs ' T1rho Slice ' int2str(sl)]; roitxt},'FontSize',16, ...
         'FontWeight','bold');
%
   if k==1
     print('-dpsc2','-r600','-fillpage',pnam2);
   else
     print('-dpsc2','-r600','-fillpage','-append',pnam2);
   end
%
   figure;
   orient landscape;
   imagesc(img2,[-100 100]);
   colormap(cmap);
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
   hb = colorbar;
   set(hb,'Limits',[10 70]);
%
   if k==1
     print('-dpsc2','-r600','-fillpage',pnam3);
   else
     print('-dpsc2','-r600','-fillpage','-append',pnam3);
   end
%
% Plot Histograms of the Cartilage T1rho Values
%
   figure;
   orient tall;
%
   hdat1 = img2(mask1);
   idv = hdat1>0&hdat1<100;            % Valid data
   hdat1 = hdat1(idv);
   n1 = size(hdat1,1);
   mn1 =  mean(hdat1);
   sd1 = std(hdat1);
   snr1 = mn1./sd1;
%
   htxt1 = ['N = ' int2str(n1)];
   htxt2 = ['Mean T1\rho = ' sprintf('%.1f',mn1)];
   htxt3 = ['SD T1\rho = ' sprintf('%.2f',sd1)];
   htxt4 = ['SNR = ' sprintf('%.3f',snr1)];
   htxt = {htxt1; htxt2; htxt3; htxt4};
%
   subplot(2,1,1);
   histogram(hdat1,25,'FaceAlpha',1,'FaceColor',[0 0 0.8]);
   axlim = axis;
   text(75,axlim(4)/2,htxt,'FontSize',11,'FontWeight','bold');
   xlabel('T1rho Values (ms)','FontSize',12,'FontWeight','bold');
   ylabel('Frequency','FontSize',12,'FontWeight','bold');
   title({[fs ' Slice ' int2str(sl) ' Histogram']; ...
          'Superficial Layer'},'FontSize',16,'FontWeight','bold');
%
   hdat2 = img2(mask2);
   idv = hdat2>0&hdat2<100;            % Valid data
   hdat2 = hdat2(idv);
   n2 = size(hdat2,1);
   mn2 =  mean(hdat2);
   sd2 = std(hdat2);
   snr2 = mn2./sd2;
%
   htxt1 = ['N = ' int2str(n2)];
   htxt2 = ['Mean T1\rho = ' sprintf('%.1f',mn2)];
   htxt3 = ['SD T1\rho = ' sprintf('%.2f',sd2)];
   htxt4 = ['SNR = ' sprintf('%.3f',snr2)];
   htxt = {htxt1; htxt2; htxt3; htxt4};
%
   subplot(2,1,2);
   histogram(hdat2,25,'FaceAlpha',1,'FaceColor',[0 0 0.8]);
   axlim = axis;
   text(75,axlim(4)/2,htxt,'FontSize',11,'FontWeight','bold');
   xlabel('T1rho Values (ms)','FontSize',12,'FontWeight','bold');
   ylabel('Frequency','FontSize',12,'FontWeight','bold');
   title({[fs ' Slice ' int2str(sl) ' Histogram']; ...
          'Deep Layer'},'FontSize',16,'FontWeight','bold');
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
% Get Cartilage T1rho Measures for Each Slice
%
% datsm contains the valid cartilage T1rho values for each bone on each
% ROI slice.
%
   for l = 1:2          % Layer
      for m = 1:3       % Bone (femur, patella or tibia)
%
         slnums(l,m,k) = sl;
         datas{l,m,k} = single(datas{l,m,k});
%
%  Cartilage T1rho Measures for Each Bone on Each Slice
%
         if ibone(k,m)&&m==1
           datsm = dats(maskf(:,l,k));
           idv = datsm>0&datsm<100;    % Valid data
           datsm = datsm(idv);
           datas{l,m,k} = datsm;
           nums(l,m,k) = size(datsm,1);
           T1rho_mn(l,m,k) = mean(datsm,1);
           T1rho_sd(l,m,k) = std(datsm,1);
           T1rho_snr(l,m,k) = T1rho_mn(l,m,k)/T1rho_sd(l,m,k);
         end
%
         if ibone(k,m)&&m==2
           datsm = dats(maskp(:,l,k));
           idv = datsm>0&datsm<100;    % Valid data
           datsm = datsm(idv);
           datas{l,m,k} = datsm;
           nums(l,m,k) = size(datsm,1);
           T1rho_mn(l,m,k) = mean(datsm,1);
           T1rho_sd(l,m,k) = std(datsm,1);
           T1rho_snr(l,m,k) = T1rho_mn(l,m,k)/T1rho_sd(l,m,k);
         end
%
         if ibone(k,m)&&m==3
           datsm = dats(maskt(:,l,k));
           idv = datsm>0&datsm<100;    % Valid data
           datsm = datsm(idv);
           datas{l,m,k} = datsm;
           nums(l,m,k) = size(datsm,1);
           T1rho_mn(l,m,k) = mean(datsm,1);
           T1rho_sd(l,m,k) = std(datsm,1);
           T1rho_snr(l,m,k) = T1rho_mn(l,m,k)/T1rho_sd(l,m,k);
         end
      end
%
% Cartilage T1rho Measures for All Bones on Each Slice
%
      slnums(l,m+1,k) = sl;            % All bones
      if l==1
        nums(l,m+1,k) = n1;
        T1rho_mn(l,m+1,k) = mn1;
        T1rho_sd(l,m+1,k) = sd1;
        T1rho_snr(l,m+1,k) = snr1;
      else
        nums(l,m+1,k) = n2;
        T1rho_mn(l,m+1,k) = mn2;
        T1rho_sd(l,m+1,k) = sd2;
        T1rho_snr(l,m+1,k) = snr2;
      end
   end
%
end
%
% Get Cartilage T1rho Measures for Each Bone Over All ROI Slices
%
data{2,3} = single([]);
na = zeros(2,3);
mna = zeros(2,3,'single');
sda = zeros(2,3,'single');
snra = zeros(2,3,'single');
%
for l = 1:2
   for m = 1:3
      dat = cell2mat(squeeze(datas(l,m,:)));
      data{l,m} = dat;
      na(l,m) = size(dat,1);
      mna(l,m) = mean(dat);
      sda(l,m) = std(dat);
      snra(l,m) = mna(l,m)/sda(l,m);
   end
end
%
% Construct Output Table for T1rho Measures on Each Slice
%
bones2 = repmat({'Femur' 'Patella' 'Tibia' 'All'},2,1);
bones2 = bones2(:);     % Bone Names
bonesa = bones2(1:6);
bones2 = repmat(bones2,nrsl,1);
%
layer = repmat((1:2)',4*nrsl,1);
lays = {'Superficial';'Deep'};
%
nums = nums(:);
idv = find(nums);       % Valid measures
%
slnums = slnums(:);
slnums = slnums(idv);
slnums = cellstr(int2str(slnums));
bones2 = bones2(idv);
layer = layer(idv);
lays2 = lays(layer);
nums = nums(idv);
T1rho_mn = T1rho_mn(:);
T1rho_mn = T1rho_mn(idv);
T1rho_sd = T1rho_sd(:);
T1rho_sd = T1rho_sd(idv);
T1rho_snr = T1rho_snr(:);
T1rho_snr = T1rho_snr(idv);
%
varnams = {'Slice#','Bone','Layer','N','Mean','SD','SNR'};
%
tab1 = table(slnums,bones2,lays2,nums,T1rho_mn,T1rho_sd,T1rho_snr, ...
           'VariableNames',varnams);
%
% Construct Output Table for T1rho Measures Over all ROI Slices
%
slnuma = cellstr(repmat('All',6,1));
laya = repmat(lays,3,1);
%
na = na(:);
mna = mna(:);
sda = sda(:);
snra = snra(:);
%
tab2 = table(slnuma,bonesa,laya,na,mna,sda,snra,'VariableNames',varnams);
%
% Combine Output Tables with T1rho Measures and Write to MS-Excel File
%
tab = [tab1; tab2];
%
writetable(tab,xnam);
%
% Save Masks, ROIS and Slice Information into MAT File
%
save(mnam,'bc','bones','f','ibone','irsl','maskb','maskf','maskp', ...
     'maskt','nrfiles','p','rois','rsl','t','tab','-append');
%
return