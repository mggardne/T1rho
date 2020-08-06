%#######################################################################
%
%                    * T1rho SLice CoMPare Program *
%
%          M-File which plots comparisons of T1rho values in the slices
%      with ROIs and plots histograms of the values over all the slices.
%
%     NOTES:  None.
%
%     05-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Load T1rho MAT Files
%
mnam = 'T1rho_voxel0.5mmslice2.0.mat';
load(mnam,'idx','irsl','mask*','series_desc','T1rhonls');
mnamp = [mnam(1:5) 'p' mnam(6:end)];
load(mnamp,'T1rhop');
%
% Get Number of Slices with ROIs and Series Description
%
nrsl = size(irsl,1);
%
fs = series_desc{idx};
fs = fs(8:end);         % Short series description w/o '3DMAPPS'
fs = strtrim(fs);       % File name description
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = zeros(128,3);    % Jet color map for cartilage
jmap(13:90,:) = jet(78);% Jet colors in scaled range 128*[0 10 70 100]./100 
jmap(1:12,:) = repmat(jmap(13,:),12,1);          % Set lower values to blue
jmap(91:128,:) = repmat(jmap(90,:),128-90,1);    % Set higher values to red
cmap = [gmap; jmap];
%
% Loop Through the Slices and Plot the T1rho Values
%
cnam = [fs '_cmp.ps'];  % Comparison print file name
%
mask = maskf|maskp|maskt;              % Combined cartilage mask
%
t1rho = cell(nrsl,2);   % T1rho for all slices (1 - UVM, 2 - Philips)
%
for k = 1:nrsl
%
   sl = irsl(k);
   msk = mask(:,k);
%
% Plot Both T1rho Values
%
   figure;
   orient tall;
   img1 = T1rhonls(:,:,sl);
   img1(img1<0) = 0;                   % Truncate data
   img1(img1>100) = 100;               % Truncate data
   img1(~msk) = img1(~msk)-100.1;
   t1rho{k,1} = img1(msk);
   subplot(2,1,1);
   imagesc(img1,[-100 100]);
   colormap(cmap);
   axis image;
   axis off;
   title([fs ' UVM T1rho Slice ' int2str(sl)],'FontSize',14, ...
         'FontWeight','bold');
   hb = colorbar;
   set(hb,'Limits',[10 70]);
%
   img2 = T1rhop(:,:,sl);
   img2(~msk) = img2(~msk)-100.1;
   t1rho{k,2} = img2(msk);
   subplot(2,1,2);
   imagesc(img2,[-100 100]);
   colormap(cmap);
   axis image;
   axis off;
   title([fs ' Philips T1rho Slice ' int2str(sl)],'FontSize',14, ...
         'FontWeight','bold');
   hb = colorbar;
   set(hb,'Limits',[10 70]);
%
   if k==1
     print('-dpsc2','-r600','-fillpage',cnam);
   else
     print('-dpsc2','-r600','-fillpage','-append',cnam);
   end
%
% Get Differences in T1rho Values
%
   msk = reshape(msk,size(img1));
   [xrng,yrng] = mask_rng(msk);
   ix = xrng(1)-2:xrng(2)+2;
   iy = yrng(1)-2:yrng(2)+2;
   mskd = msk(ix,iy);
   imgd = img1(ix,iy)-img2(ix,iy);
   imgd(~mskd) = NaN;
%
% Plot Differences in T1rho Values
%
   figure;
   orient tall;
%
   subplot(2,1,1);
   imagesc(imgd,[-0.01 0.025]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' T1rho Differences Slice ' int2str(sl)],'FontSize',14, ...
         'FontWeight','bold');
   colorbar;
%
   x = linspace(-1e-03,2.5e-2,25)';
%
   subplot(2,1,2);
   h = histogram(imgd(:),x,'FaceAlpha',1,'FaceColor',[0 0 0.8]);
%
   xlabel('T1rho Differences (ms)','FontSize',11,'FontWeight','bold');
   ylabel('Frequency','FontSize',11,'FontWeight','bold');
   title([fs ' Slice ' int2str(sl) ' Histogram'],'FontSize',14, ...
         'FontWeight','bold');
%
   print('-dpsc2','-r600','-fillpage','-append',cnam);

%
end
%
% Get T1rho Values for All the Slices
%
t1ru = cell2mat(t1rho(:,1));           % UVM values
t1rp = cell2mat(t1rho(:,2));           % Philips values
t1ru(t1ru<0) = NaN;     % Values outside of cartilage
t1rp(t1rp<0) = NaN;     % Values outside of cartilage
t1rd = t1ru-t1rp;       % Differences in cartilage T1rho values
%
% Plot Overall Cartilage Histograms
%
figure;
orient tall;
%
subplot(3,1,1)
histogram(t1ru,20,'FaceColor',[0 0 0.8],'FaceAlpha',1);
%
xlabel('T1rho (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title([fs ' Combined UVM T1rho Histogram'],'FontSize',12, ...
      'FontWeight','bold');
%
subplot(3,1,2)
histogram(t1rp,20,'FaceColor',[0 0 0.8],'FaceAlpha',1);
%
xlabel('T1rho (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title([fs ' Combined Philips T1rho Histogram'],'FontSize',12, ...
      'FontWeight','bold');
%
subplot(3,1,3)
histogram(t1rd,x,'FaceColor',[0 0 0.8],'FaceAlpha',1);
%
xlabel('T1rho Differences (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title([fs ' Combined T1rho Differences Histogram'],'FontSize',12, ...
      'FontWeight','bold');
%
print('-dpsc2','-r600','-fillpage','-append',cnam);
%
return