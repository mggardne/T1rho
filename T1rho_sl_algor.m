%#######################################################################
%
%                   * T1rho SLice ALGORithm Program *
%
%          M-File which plots results about the algorithm for computing
%      T1rho for different regions.
%
%     NOTES:  None.
%
%     07-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Load T1rho MAT Files
%
mnam = 'T1rho_voxel0.5mmslice2.0.mat';
load(mnam,'exit_flag','fnams','id','idx','irsl','mask*','nls_amp', ...
          'nls_nres','npx','nslt','rsl','series_desc','T1rhonls');
%
dir_nam = fullfile('v3','DICOM');
%
fnams = fnams(id);      % File names for this series
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
gmap = gray(128);       % Gray color map for pixels not within mask
jmap = jet(128);        % Jet color map for pixels within mask
cmap = [gmap; jmap];
%
% Set Up Variables and Arrays for Loop
%
snam = [fs '_algor.ps']; % Algorithm print file name
hnam = [fs '_hist.ps'];  % Histogram print file name
%
maskc = maskf|maskp|maskt;             % Combined cartilage mask (femur, patella, tibia)
%
na = sum(maska);        % Number of pixels in mask per slice
nat = sum(na);          % Total number of pixels in mask for all ROI slices
nb = sum(maskb);        % Number of pixels in mask per slice
nbt = sum(nb);          % Total number of pixels in mask for all ROI slices
nc = sum(maskc);        % Number of pixels in mask per slice
nct = sum(nc);          % Total number of pixels in mask for all ROI slices
%
rslm = rsl+repmat(0:nslt-1,4,1);       % Slices with all the spin lock times
%
t1rhoa = zeros(nat,1,'single');        % T1rho for air for all slices
t1rhob = zeros(nbt,1,'single');        % T1rho for bone for all slices
t1rhoc = zeros(nct,1,'single');        % T1rho for cartilage for all slices
%
nls_ampa = zeros(nat,1,'single');      % Fit amplitude
nls_ampb = zeros(nbt,1,'single');
nls_ampc = zeros(nct,1,'single');
%
nls_nresa = zeros(nat,1,'single');     % Fit normalized residuals
nls_nresb = zeros(nbt,1,'single');
nls_nresc = zeros(nct,1,'single');
%
exit_flaga = zeros(nat,1,'uint8');     % Fit exit flag
exit_flagb = zeros(nbt,1,'uint8');
exit_flagc = zeros(nct,1,'uint8');
%
spin0a = zeros(nat,1,'single');        % T1 signal at spin lock time = 0
spin0b = zeros(nbt,1,'single');
spin0c = zeros(nct,1,'single');
%
sse = zeros(npx,npx,nslt);
sse1 = zeros(npx,npx);
%
sse_spina = zeros(nat,1,'single');     % Sum of squared errors
sse_spinb = zeros(nbt,1,'single');
sse_spinc = zeros(nct,1,'single');
%
% Loop Through the Slices and Plot T1rho Values
%
for k = 1:nrsl

%
%    slk = rsl(k);        % Slice number in spin lock images
   sl = irsl(k);        % Slice number in T1rho images
   sll = int2str(sl);
%
   ida = sum(na(1:k));
   ida = ida-na(k)+1:ida;
   idb = sum(nb(1:k));
   idb = idb-nb(k)+1:idb;
   idc = sum(nc(1:k));
   idc = idc-nc(k)+1:idc;
%
   ma = maska(:,k);
   mb = maskb(:,k);
   mc = maskc(:,k);
%
% Get the Four Outcomes From the Three Areas
%
   t1 = T1rhonls(:,:,sl);
   t1rhoa(ida) = t1(ma);
   t1rhob(idb) = t1(mb);
   t1rhoc(idc) = t1(mc);
%
   a1 = nls_amp(:,:,sl);
   nls_ampa(ida) = a1(ma);
   nls_ampb(idb) = a1(mb);
   nls_ampc(idc) = a1(mc);
%
   r1 = nls_nres(:,:,sl);
   nls_nresa(ida) = r1(ma);
   nls_nresb(idb) = r1(mb);
   nls_nresc(idc) = r1(mc);
%
   ef1 = exit_flag(:,:,sl);
   exit_flaga(ida) = ef1(ma);
   exit_flagb(idb) = ef1(mb);
   exit_flagc(idc) = ef1(mc);
%   
   figure;
   orient tall;
   subplot(3,2,1);
   imagesc(t1,[-10 200]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' T1\rho Slice ' sll],'FontSize',12, ...
         'FontWeight','bold');
   hb = colorbar;
%   set(hb,'Limits',[10 70]);
%
   subplot(3,2,2);
   imagesc(a1,[0 200]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' Amplitude Slice ' sll],'FontSize',12, ...
         'FontWeight','bold');
   hb = colorbar;
%   set(hb,'Limits',[0 200]);
%
   subplot(3,2,3);
   imagesc(r1,[0 1]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' Residuals Slice ' sll],'FontSize',12, ...
         'FontWeight','bold');
   hb = colorbar;
%   set(hb,'Limits',[0 1]);
%
   subplot(3,2,4);
   imagesc(ef1,[-2 4]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' Exit Flags Slice ' sll],'FontSize',12, ...
         'FontWeight','bold');
   hb = colorbar;
   set(hb,'Limits',[-2 4]);
%
% Initialize Arrays
%
   ssea = zeros(na(k),nslt,'single');
   sseb = zeros(nb(k),nslt,'single');
   ssec = zeros(nc(k),nslt,'single');
%
% Loop through Spin Lock Times
%
   for l = 1:nslt       % Loop through spin lock times
%
      n = rslm(k,l);
      fnam = fnams{n};  % Filename for this spin lock time
      fnam = fullfile(dir_nam,fnam);
      fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                 ', Slice:  ' sll ', Spin lock time # ' int2str(l)]);
%
% Load and Scale Image
%
      img = dicomread(fnam);
      info = dicominfo(fnam);
      sl = single(info.RescaleSlope);
      offst = single(info.RescaleIntercept);     % Usually zero
      img = single(img);
      r = (img-offst)/sl;
      ra = r(ma);
      rb = r(mb);
      rc = r(mc);
      if l==1
        spin0a(ida) = ra;
        spin0b(idb) = rb;
        spin0c(idc) = rc;
%
        subplot(3,2,5)
        imagesc(r,[0 200]);
        colormap jet;
        axis image;
        axis off;
        title([fs ' T1 Spin Lock 0 ms Slice ' sll],'FontSize',12, ...
              'FontWeight','bold');
        hb = colorbar;
        set(hb,'Limits',[0 200]);
%
      end
%
      sse(:,:,l) = r;
      ssea(:,l) = ra;
      sseb(:,l) = rb;
      ssec(:,l) = rc;
%
   end
%
   fprintf(1,'\n');     % Extra line between sets of spin lock times
%
% Calculate Sum of Squares of T1 Values
%
   mn_ssea = mean(ssea,2);
   ssea = ssea-repmat(mn_ssea,1,nslt);
   sse_spina(ida) = sum(ssea.*ssea,2);
%
   mn_sseb = mean(sseb,2);
   sseb = sseb-repmat(mn_sseb,1,nslt);
   sse_spinb(idb) = sum(sseb.*sseb,2);
%
   mn_ssec = mean(ssec,2);
   ssec = ssec-repmat(mn_ssec,1,nslt);
   sse_spinc(idc) = sum(ssec.*ssec,2);
%
% Plot SSE of T1 Values
%
   mn_sse = mean(sse,3);
   sse = sse-repmat(mn_sse,1,1,nslt);
   sse1 = sum(sse.*sse,3);
%
   subplot(3,2,6);
   imagesc(sse1,[0 25000]);
   colormap jet;
   axis image;
   axis off;
   title([fs ' SSE Slice ' sll],'FontSize',12, ...
         'FontWeight','bold');
   hb = colorbar;
%   set(hb,'Limits',[0 25000]);
%
   if k==1
     print('-dpsc2','-r600','-fillpage',snam);
   else
     print('-dpsc2','-r600','-fillpage','-append',snam);
   end
%
end
%
% Get Summary Statistics for All the Slices
%
fprintf(1,'\nRange of T1rho in air:  %g-%g\n',min(t1rhoa),max(t1rhoa));
fprintf(1,'Mean of T1rho in air:  %g\n',mean(t1rhoa));
fprintf(1,'Median of T1rho in air:  %g\n',median(t1rhoa));
fprintf(1,'SD of T1rho in air:  %g\n\n',std(t1rhoa));
fprintf(1,'Range of T1rho in bone:  %g-%g\n',min(t1rhob),max(t1rhob));
fprintf(1,'Mean of T1rho in bone:  %g\n',mean(t1rhob));
fprintf(1,'Median of T1rho in bone:  %g\n',median(t1rhob));
fprintf(1,'SD of T1rho in bone:  %g\n\n',std(t1rhob));
fprintf(1,'Range of T1rho in cartilage:  %g-%g\n',min(t1rhoc), ...
        max(t1rhoc));
fprintf(1,'Mean of T1rho in cartilage:  %g\n',mean(t1rhoc));
fprintf(1,'Median of T1rho in cartilage:  %g\n',median(t1rhoc));
fprintf(1,'SD of T1rho in cartilage:  %g\n\n',std(t1rhoc));
%
fprintf(1,'Range of Residuals in air:  %g-%g\n',min(nls_nresa), ...
        max(nls_nresa));
fprintf(1,'Mean of Residuals in air:  %g\n',mean(nls_nresa));
fprintf(1,'Median of Residuals in air:  %g\n',median(nls_nresa));
fprintf(1,'SD of Residuals in air:  %g\n\n',std(nls_nresa));
fprintf(1,'Range of Residuals in bone:  %g-%g\n',min(nls_nresb), ...
        max(nls_nresb));
fprintf(1,'Mean of Residuals in bone:  %g\n',mean(nls_nresb));
fprintf(1,'Median of Residuals in bone:  %g\n',median(nls_nresb));
fprintf(1,'SD of Residuals in bone:  %g\n\n',std(nls_nresb));
fprintf(1,'Range of Residuals in cartilage:  %g-%g\n', ...
        min(nls_nresc),max(nls_nresc));
fprintf(1,'Mean of Residuals in cartilage:  %g\n',mean(nls_nresc));
fprintf(1,'Median of Residuals in cartilage:  %g\n',median(nls_nresc));
fprintf(1,'SD of Residuals in cartilage:  %g\n\n',std(nls_nresc));
%
fprintf(1,'Range of Amplitude in air:  %g-%g\n',min(nls_ampa), ...
        max(nls_ampa));
fprintf(1,'Mean of Amplitude in air:  %g\n',mean(nls_ampa));
fprintf(1,'Median of Amplitude in air:  %g\n',median(nls_ampa));
fprintf(1,'SD of Amplitude in air:  %g\n\n',std(nls_ampa));
fprintf(1,'Range of Amplitude in bone:  %g-%g\n',min(nls_ampb), ...
        max(nls_ampb));
fprintf(1,'Mean of Amplitude in bone:  %g\n',mean(nls_ampb));
fprintf(1,'Median of Amplitude in bone:  %g\n',median(nls_ampb));
fprintf(1,'SD of Amplitude in bone:  %g\n\n',std(nls_ampb));
fprintf(1,'Range of Amplitude in cartilage:  %g-%g\n', ...
        min(nls_ampc),max(nls_ampc));
fprintf(1,'Mean of Amplitude in cartilage:  %g\n',mean(nls_ampc));
fprintf(1,'Median of Amplitude in cartilage:  %g\n',median(nls_ampc));
fprintf(1,'SD of Amplitude in cartilage:  %g\n\n',std(nls_ampc));
%
fprintf(1,'Range of T1 Spin Lock 0 ms in air:  %g-%g\n',min(spin0a), ...
        max(spin0a));
fprintf(1,'Mean of T1 Spin Lock 0 ms in air:  %g\n',mean(spin0a));
fprintf(1,'Median of T1 Spin Lock 0 ms in air:  %g\n',median(spin0a));
fprintf(1,'SD of T1 Spin Lock 0 ms in air:  %g\n\n',std(spin0a));
fprintf(1,'Range of T1 Spin Lock 0 ms in bone:  %g-%g\n', ...
        min(spin0b),max(spin0b));
fprintf(1,'Mean of T1 Spin Lock 0 ms in bone:  %g\n',mean(spin0b));
fprintf(1,'Median of T1 Spin Lock 0 ms in bone:  %g\n',median(spin0b));
fprintf(1,'SD of T1 Spin Lock 0 ms in bone:  %g\n\n',std(spin0b));
fprintf(1,'Range of T1 Spin Lock 0 ms in cartilage:  %g-%g\n', ...
        min(spin0c),max(spin0c));
fprintf(1,'Mean of T1 Spin Lock 0 ms in cartilage:  %g\n',mean(spin0c));
fprintf(1,'Median of T1 Spin Lock 0 ms in cartilage:  %g\n', ...
        median(spin0c));
fprintf(1,'SD of T1 Spin Lock 0 ms in cartilage:  %g\n\n',std(spin0c));
%
fprintf(1,'Range of SSE in air:  %g-%g\n',min(sse_spina), ...
        max(sse_spina));
fprintf(1,'Mean of SSE in air:  %g\n',mean(sse_spina));
fprintf(1,'Median of SSE in air:  %g\n',median(sse_spina));
fprintf(1,'SD of SSE in air:  %g\n\n',std(sse_spina));
fprintf(1,'Range of SSE in bone:  %g-%g\n',min(sse_spinb), ...
        max(sse_spinb));
fprintf(1,'Mean of SSE in bone:  %g\n',mean(sse_spinb));
fprintf(1,'Median of SSE in bone:  %g\n',median(sse_spinb));
fprintf(1,'SD of SSE in bone:  %g\n\n',std(sse_spinb));
fprintf(1,'Range of SSE in cartilage:  %g-%g\n',min(sse_spinc), ...
        max(sse_spinc));
fprintf(1,'Mean of SSE in cartilage:  %g\n',mean(sse_spinc));
fprintf(1,'Median of SSE in cartilage:  %g\n',median(sse_spinc));
fprintf(1,'SD of SSE in cartilage:  %g\n\n',std(sse_spinc));
%
fprintf(1,'Range of Exit Flags in air:  %g-%g\n',min(exit_flaga), ...
        max(exit_flaga));
fprintf(1,'Mean of Exit Flags in air:  %g\n',mean(exit_flaga));
fprintf(1,'Median of Exit Flags in air:  %g\n',median(exit_flaga));
fprintf(1,'SD of Exit Flags in air:  %g\n\n',std(single(exit_flaga)));
fprintf(1,'Range of Exit Flags in bone:  %g-%g\n',min(exit_flagb), ...
        max(exit_flagb));
fprintf(1,'Mean of Exit Flags in bone:  %g\n',mean(exit_flagb));
fprintf(1,'Median of Exit Flags in bone:  %g\n',median(exit_flagb));
fprintf(1,'SD of Exit Flags in bone:  %g\n\n',std(single(exit_flagb)));
fprintf(1,'Range of Exit Flags in cartilage:  %g-%g\n', ...
        min(exit_flagc),max(exit_flagc));
fprintf(1,'Mean of Exit Flags in cartilage:  %g\n',mean(exit_flagc));
fprintf(1,'Median of Exit Flags in cartilage:  %g\n', ...
        median(exit_flagc));
fprintf(1,'SD of Exit Flags in cartilage:  %g\n\n', ...
        std(single(exit_flagc)));
%
% Plot Overall Histograms
%
figure;
orient tall;
%
subplot(3,2,1);
histogram(t1rhoa,[-0.15 -0.05 0.05 0.15],'FaceColor', ...
          [0 0 0.8],'FaceAlpha',1);
%
xlabel('T1\rho in Air (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1\rho in Air','FontSize',12,'FontWeight','bold');
%
subplot(3,2,3);
histogram(t1rhob,[-Inf -200:100:2000 Inf],'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('T1\rho in Bone (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1\rho in Bone','FontSize',12,'FontWeight','bold');
%
subplot(3,2,5);
histogram(t1rhoc,[-Inf -20:10:200 +Inf],'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('T1\rho in Cartilage (ms)','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1\rho in Cartilage','FontSize',12,'FontWeight','bold');
%
subplot(3,2,2);
histogram(nls_nresa,[-0.15 -0.05 0.05 0.15],'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Residuals in Air','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Normalized Residuals in Air','FontSize',12,'FontWeight','bold');
%
subplot(3,2,4);
histogram(nls_nresb,0:0.05:0.8,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Residuals in Bone','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Normalized Residuals in Bone','FontSize',12,'FontWeight','bold');
%
subplot(3,2,6);
histogram(nls_nresc,0:0.01:0.8,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Residuals in Cartilage','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Normalized Residuals in Cartilage','FontSize',12, ...
      'FontWeight','bold');
%
print('-dpsc2','-r600','-fillpage',hnam);
%
figure;
orient tall;
%
subplot(3,2,1);
histogram(nls_ampa,[-0.15 -0.05 0.05 0.15],'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Amplitude in Air','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Amplitude in Air','FontSize',12,'FontWeight','bold');
%
subplot(3,2,3);
histogram(nls_ampb,0:2.5:125,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Amplitude in Bone','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Amplitude in Bone','FontSize',12,'FontWeight','bold');
%
subplot(3,2,5);
histogram(nls_ampc,0:5:200,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Amplitude in Cartilage','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Amplitude in Cartilage','FontSize',12,'FontWeight','bold');
%
subplot(3,2,2);
histogram(spin0a,0:0.1:0.9,'FaceColor', ...
          [0 0 0.8],'FaceAlpha',1);
%
xlabel('T1 Spin Lock 0 ms in Air','FontSize',10,'FontWeight', ...
       'bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1 Spin Lock 0 ms in Air','FontSize',12,'FontWeight','bold');
%
subplot(3,2,4);
histogram(spin0b,0:2.5:125,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('T1 Spin Lock 0 ms in Bone','FontSize',10,'FontWeight', ...
       'bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1 Spin Lock 0 ms in Bone','FontSize',12,'FontWeight','bold');
%
subplot(3,2,6);
histogram(spin0c,0:5:200,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('T1 Spin Lock 0 ms in Cartilage','FontSize',10, ...
       'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('T1 Spin Lock 0 ms in Cartilage','FontSize',12,'FontWeight', ...
      'bold');
%
print('-dpsc2','-r600','-fillpage','-append',hnam);
%
figure;
orient tall;
%
subplot(3,2,1);
histogram(sse_spina,0:0.05:0.7,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('SSE in Air','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Sum of Squared Errors in Air','FontSize',12, ...
      'FontWeight','bold');
%
subplot(3,2,3);
histogram(sse_spinb,0:500:12000,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('SSE in Bone','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Sum of Squared Errors in Bone','FontSize',12, ...
      'FontWeight','bold');
%
subplot(3,2,5);
histogram(sse_spinc,0:500:24000,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('SSE in Cartilage','FontSize',10,'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Sum of Squared Errors in Cartilage','FontSize',12, ...
      'FontWeight','bold');
%
subplot(3,2,2);
histogram(exit_flaga,[-0.15 -0.05 0.05 0.15],'FaceColor', ...
          [0 0 0.8],'FaceAlpha',1);
%
xlabel('Exit Flags in Air','FontSize',10,'FontWeight', ...
       'bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Exit Flags in Air','FontSize',12,'FontWeight','bold');
%
subplot(3,2,4);
histogram(exit_flagb,-0.5:1:4.5,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Exit Flags in Bone','FontSize',10,'FontWeight', ...
       'bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Exit Flags in Bone','FontSize',12,'FontWeight','bold');
%
subplot(3,2,6);
histogram(exit_flagc,-0.5:1:4.5,'FaceColor',[0 0 0.8], ...
          'FaceAlpha',1);
%
xlabel('Exit Flags in Cartilage','FontSize',10, ...
       'FontWeight','bold');
ylabel('Frequency','FontSize',10,'FontWeight','bold');
title('Exit Flags in Cartilage','FontSize',12,'FontWeight', ...
      'bold');
%
print('-dpsc2','-r600','-fillpage','-append',hnam);
%
return