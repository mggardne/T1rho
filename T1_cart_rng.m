%#######################################################################
%
%                    * T1 CARTilage RaNGe Program *
%
%          M-File which reads knee MRI data with different spin lock
%     times and plots the data for the cartilage regions.  The 
%     superficial and deep layers are plotted separately.  Also
%     generates a box and whiskers plot of the combined cartilage data.
%
%     NOTES:  None.
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
dat1 = cell(1,nrsl);    % T1 values in layer 1
dat2 = cell(1,nrsl);    % T1 values in layer 2
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
   id1 = find(mask1(:,k));             % Index to layer 1 (superficial)
   idc1 = ismember(idp,id1);
   np(k) = size(idp,1); % Number of pixels on this slice
   npk = np(k);
%
% Initialize Arrays for Loop
%
   dat2d = zeros(nslt,npk);            % Data for all spin lock times
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
   dat1{k} = dat2d(:,idc1);
   dat2{k} = dat2d(:,~idc1);
%
end
%
dat1 = cell2mat(dat1)'; % Combine slice data
dat2 = cell2mat(dat2)'; % Combine slice data
dat = [dat1; dat2];     % Combine layer data
%
% Plot Data for Each Spin Lock Time
%
figure;
orient landscape;
h1 = plot(slt'-0.5,dat1,'b.');
hold on;
h2 = plot(slt'+0.5,dat2,'r.');
set(gca,'XTick',[0 5 10 20:20:80]);
xlabel('Spin Lock Time (ms)','FontSize',12,'FontWeight','bold');
ylabel('T1','FontSize',12,'FontWeight','bold');
title([fs ' Cartilage'],'FontSize',16,'FontWeight','bold');
axlim = axis;
plot([-5 85],[0 0],'k:');
axis([-5 85 -5 axlim(4)]);
legend([h1(1); h2(1)],{'Superficial'; 'Deep'},'FontSize',11, ...
       'Location','northeast');
%
pnam = 'T1_cart_rng.ps';
print('-dpsc2','-r600','-fillpage',pnam);
%
% Box and Whiskers Plot
%
figure;
orient landscape;
h = boxplot(dat,'positions',slt','labels',cellstr(int2str(slt)));
hold on;
xlabel('Spin Lock Time (ms)','FontSize',12,'FontWeight','bold');
ylabel('T1','FontSize',12,'FontWeight','bold');
title({[fs ' Cartilage']; 'Box and Whiskers Plot'},'FontSize',16, ...
      'FontWeight','bold');
axlim = axis;
plot(axlim(1:2),[0 0],'k:');
%
print('-dpsc2','-r600','-fillpage','-append',pnam);
%
return