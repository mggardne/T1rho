%#######################################################################
%
%                      * T1 CARTilage SD Program *
%
%          M-File which reads knee MRI data with different spin lock
%     times and plots the standard deviations for the cartilage regions.
%     The different spin lock times are plotted separately.
%
%     NOTES:  None.
%
%     08-Sep-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Pick MRI Series to Analyze
%
mnams = dir('T1rho_*.mat');
mnams = {mnams.name}';
%
idm = menu('Pick a MAT File to Analyze',mnams);
mnam = mnams{idm};
%
fs = extractAfter(mnam,'_');
fs = fs(1:end-4);
%
% Get DICOM File Names, Masks, Spin Lock Times and Slice Numbers
%
load(mnam,'fnams','irsl','maskf','maskp','maskt','npx','nslt', ...
     'rsl','slt');
nrsl = size(rsl,1);     % Number of slices with cartilage ROIs
npx2 = npx*npx;         % Number of elements in an image
%
% Combine Masks to Get Mask for All Cartilage on a Slice
%
mask1 = squeeze(maskf(:,1,:)|maskp(:,1,:)|maskt(:,1,:));   % Layer 1 (superficial)
mask2 = squeeze(maskf(:,2,:)|maskp(:,2,:)|maskt(:,2,:));   % Layer 2 (deep)
maskc = mask1|mask2;    % All cartilage mask
%
% Setup Color Map and PS File Names
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
pnam1 = ['T1_SD_' fs '.ps'];           % SD PS print file name
pnam2 = ['T1_SNR_' fs '.ps'];          % SNR PS print file name
%
SDtxt = [fs ' Standard Deviations'];
SNRtxt = [fs ' Signal to Noise Ratios'];
%
% Loop through Slices
%
for k = 1:nrsl
%
% Slice Information
%
   n = rsl(k);          % Slice number
%
% Get Mask for this Slice
%
   mask = maskc(:,k);
   id = find(mask);
   npts = size(id,1);
%
% Initialize Arrays for Spin Lock Times Loop
%
   dat2d = zeros(npx*npx,nslt);        % Data for all spin lock times
   rmn = zeros(npts,nslt);             % Means for all spin lock times
   sd = zeros(npts,nslt);              % SD for all spin lock times
   snr = zeros(npts,nslt);             % Signal-to-noise ratio (SNR) for all spin lock times
%
% Loop through Spin Lock Times
%
   for l = 1:nslt       % Loop through spin lock times
%
      nf = n-1+l;
      sll = int2str(nf);               % Slice number as letters
      fnam = fnams{nf}; % Filename for this spin lock time
      sltl = int2str(slt(l));
      fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                 ', Slice:  ' sll ', Spin lock time:  ' sltl ' ms']);
%
% Load and Scale Slice Image
%
      img = dicomread(fnam);
      img = single(img);
      info = dicominfo(fnam);
      sl = single(info.RescaleSlope);
      offst = single(info.RescaleIntercept);     % Usually zero
      img = double((img-offst)./sl);
      dat2d(:,l) = img(:);
      mx = max(img(:))+0.1;
      img = (img-mx)./mx;
%
% Calculate Means, Standard Deviations and SNR within Cartilage
%
      [idnn,idv] = get_nn_idx(id,npx,true);
      [rmn(:,l),sd(:,l)] = stat_nan(dat2d(:,l),idnn,idv);
      snr(:,l) = rmn(:,l)./sd(:,l);
%
% Get Range for Approximately 98% of Values
%
      cutoff = 0.98;
      [nsd,edg] = histcounts(sd(:,l),200);
      nsd = cumsum(nsd)./npts;
      idsd = find(nsd>=cutoff,1,'first');
      mxsd = edg(idsd); % Lower edge of bin
      mxsd = ceil(mxsd);
%
      [nsnr,edg] = histcounts(snr(:,l),200);
      nsnr = cumsum(nsnr)./npts;
      idsnr = find(nsnr>=cutoff,1,'first');
      mxsnr = edg(idsnr);              % Lower edge of bin
      mxsnr = 10*ceil(mxsnr/10);
%
% Plot Slice SD at this Spin Lock Time
%
      img1 = img;
      img1(id) = sd(:,l)./mxsd;
%
      figure;
      orient landscape;
      imagesc(img1,[-1 1]);
      colormap(cmap);
      axis image;
      axis off;
      title({SDtxt; ['T1 Slice ' sll ', ' sltl, ...
            ' ms Spin Lock Time']},'FontSize',16,'FontWeight','bold');
      hb = colorbar;
      set(hb,'Limits',[0 1]);
      ticks = hb.Ticks';
      ticks = ceil(mxsd)*ticks;
      hb.TickLabels = cellstr(num2str(ticks,'%g'));
%
% Print SD Plot to PS File at this Spin Lock Time
%
      if k==1&&l==1
        print('-dpsc2','-r600','-fillpage',pnam1);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam1);
      end
%
% Plot Slice SNR at this Spin Lock Time
%
      img2 = img;
      img2(id) = snr(:,l)./mxsnr;
%
      figure;
      orient landscape;
      imagesc(img2,[-1 1]);
      colormap(cmap);
      axis image;
      axis off;
      title({SNRtxt; ['T1 Slice ' sll ', ' sltl, ... 
            ' ms Spin Lock Time']},'FontSize',16,'FontWeight','bold');
      hb = colorbar;
      set(hb,'Limits',[0 1]);
      ticks = hb.Ticks';
      ticks = 10*ceil(mxsnr/10)*ticks;
      hb.TickLabels = cellstr(num2str(ticks,'%g'));
%
% Print SNR Plot to PS File at this Spin Lock Time
%
      if k==1&&l==1
        print('-dpsc2','-r600','-fillpage',pnam2);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam2);
      end
%
   end
%
   fprintf(1,'\n');     % Line between slices
%
end
%
return