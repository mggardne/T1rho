%#######################################################################
%
%                    * T1 CARTilage RaNGe Program *
%
%          M-File which reads knee MRI data with different spinlock
%     times and plots the data for the cartilage regions.  


%     fits.  The T1rho maps and square root of the sum of squared
%     errors are saved to Matlab MAT files.  Slice and overall
%     statistics are written to a MS-Excel spreadsheet.
%
%     NOTES:  1.  Must have M-files exp_fun1.m and ncombi.m in the
%             current directory or path.
%
%     04-Sep-2020 * Mack Gardner-Morse
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
load(mnam,'fnams','irsl','maskf','maskp','maskt','nslt','rsl','slt');
nrsl = size(rsl,1);     % Number of slices with cartilage ROIs
%
% Combine Masks to Get Mask for All Cartilage on a Slice
%
mask1 = squeeze(maskf(:,1,:)|maskp(:,1,:)|maskt(:,1,:));   % Layer 1 (superficial)
mask2 = squeeze(maskf(:,2,:)|maskp(:,2,:)|maskt(:,2,:));   % Layer 2 (deep)
maskc = mask1|mask2;    % All cartilage mask
%
% Set Up Arrays for Loop
%
np = zeros(nrsl,1);     % Number of cartilage pixels per slice
%
dat = cell(nrsl,1);     % T1 values
%
% Loop through Slices
%
for k = 1:nrsl
%
% Slice Information
%
   n = rsl(k);          % Slice number
   sll = int2str(n);    % Slice number as letters
%
% Cartilage Mask
%  
   mask = maskc(:,k);
   idp = find(mask);    % Index to cartilage pixels on this slice
   np(k) = size(idp,1); % Number of pixels on this slice
   npk = np(k);
%
% Initialize Arrays for Loop
%
   dat2d = zeros(nslt,npk);            % Data for all spin lock times
%
% Twelve Random Pixels to Plot
%
   ipp = unique(floor(rand(npk,1)*npk+1),'stable');
   ipp = sort(ipp(1:12));              % Index to pixels to plot
   ipps = int2str(idp(ipp));           % Pixel number as a string
%
% Loop through Spin Lock Times
%
   for l = 1:nslt       % Loop through spin lock times
%
      nf = n-1+l;
      fnam = fnams{nf}; % Filename for this spin lock time
      fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                 ', Slice:  ' sll ', Spin lock time:  ' ...
                 int2str(slt(l)) ' ms']);
%
% Load and Scale Slice Image
%
      img = dicomread(fnam);
      img = img(mask);  % Get just cartilage pixels
      info = dicominfo(fnam);
      sl = single(info.RescaleSlope);
      offst = single(info.RescaleIntercept);     % Usually zero
      img = single(img');
      dat2d(l,:) = double((img-offst)/sl);
   end
%
   fprintf(1,'\n');     % Line between slices
%
   dat{k} = dat2d;
%
end
%
dat = cell2mat(dat')';  % Combine slice data
%
% Plot Data for Each Spin Lock Time
%
figure;
orient landscape;
plot(slt',dat,'k.');
hold on;
xlabel('Spin Lock Time (ms)','FontSize',12,'FontWeight','bold');
ylabel('T1','FontSize',12,'FontWeight','bold');
title([fs ' Cartilage'],'FontSize',16,'FontWeight','bold');
axlim = axis;
plot([-5 85],[0 0],'k:');
axis([-5 85 -5 axlim(4)]);
%
pnam = 'T1_cart_rng.ps';
print('-dpsc2','-r600','-fillpage',pnam);
%
% Box and Whiskers Plot
%
figure;
orient landscape;
h = boxplot(dat,'positions',slt','labels',cellstr(int2str(slt)));
xlabel('Spin Lock Time (ms)','FontSize',12,'FontWeight','bold');
ylabel('T1','FontSize',12,'FontWeight','bold');
title({[fs ' Cartilage']; 'Box and Whiskers Plot'},'FontSize',16, ...
      'FontWeight','bold');
%
print('-dpsc2','-r600','-fillpage','-append',pnam);
%
return