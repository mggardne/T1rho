%#######################################################################
%
%                       * T1rho_map 2 Program *
%
%          M-File which reads knee MRI data with different spinlock
%     times and fits an exponential to the data using nonlinear least
%     squares.  The T1rho map and normalized residuals are output to
%     DICOM format files.
%
%     NOTES:  1.  The MRI data files must be in the subdirectory
%             v3\DICOM\.
%
%             2.  Must have M-File exp_fun1.m in the current directory
%             or path.
%
%     10-Jul-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Print Parameter
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots
%
% Nonlinear Least Squares Parameter
%
nskip = false;          % Don't skip nonlinear calculations
% nskip = true;           % Skip nonlinear calculations
%
% Number of Spin Lock Times and Pixels
%
nslt = 7;               % Number of spin lock times
slt = [0 5 10 20:20:80]';              % Spin lock times
npx = 512;              % Number of pixels in each direction (symmetrical)
npx2 = npx*npx;         % Number of pixels in slice
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on');
%
T1rho0 = 80;            % Initial value for T1rho
fun = @exp_fun1;        % Exponential function
%
% Get Data Subdirectory and Filenames for the Four Series
%
dir_nam = fullfile('v3','DICOM');
fnams = dir(fullfile(dir_nam,'IM_*'));
fnams = {fnams.name}';
nfiles = size(fnams,1);
%
% Get Index to the Four Series Based on the Numbers in the File Names
%
a = [ ['['; repmat(' ',nfiles-1,1)] char(extractAfter(fnams,'_')) ...
    [repmat(';',nfiles-1,1); ']'] ];
a = a';
a = a(:)';
fnums = eval(a);        % Numbers in the filenames
fdif = diff(fnums);
id1 = find(fdif>1);
id1 = [id1(1); id1(2:end); nfiles];    % Index to series
%
% Get Series Descriptions
%
series_desc = cell(4,1);
for k = 1:4
   id = id1(k)+1;
   info = dicominfo(fullfile(dir_nam,fnams{id}));
   series_desc{k} = info.SeriesDescription;
end
%
% Pick Series and Get Index to Files
%
idx = menu('Pick a Series to Analyze',series_desc);
id = id1(idx)+1:id1(idx+1);            % Index to filenames
nf = length(id);        % Number of files
nsl = nf/nslt;          % Number of slices
plt_sl = floor(nsl/4);
plt_sl = plt_sl:plt_sl:3*plt_sl;       % Slices to plot
%
%  Get File Name Description and DICOM Header Information
%
fs = series_desc{idx};
fs = fs(8:end);         % Short series description w/o '3DMAPPS'
fs = strtrim(fs);       % File name description
ds = ['T1rhoS' int2str(idx+2) ' ' fs]; % New DICOM description with Series number
maxint = 2^16-1;        % Maximum unit16
rslope = 100/maxint;    % Slope to scale T1rho values
ressl = 1/maxint;       % Slope for normalized residual values
if iprt
  pnam = ['T1rho_' fs '.ps'];          % Plot file name
end

%
% Set Up Arrays for Loop
%
xdat = [slt ones(nslt,1)];             % Spin lock times in ms and exponential amplitude
dat3d = zeros(npx,npx,nslt);           % 2D slice array with all the spin lock times
zmat = false(npx,npx,nsl);             % Pixels with small (zero) values
%
T1rholls = zeros(npx,npx,nsl,'single');% Linear least squares time constant
lls_amp = zeros(npx,npx,nsl,'single'); % Linear least squares amplitudes
lls_nres = zeros(npx,npx,nsl,'single');% Linear least squares normalized residuals
%
T1rhonls = zeros(npx,npx,nsl,'single');% Nonlinear least squares time constant
nls_amp = zeros(npx,npx,nsl,'single'); % Nonlinear least squares amplitude
nls_nres = zeros(npx,npx,nsl,'single');% Nonlinear least squares normalized residuals
exit_flag = zeros(npx,npx,nsl,'int8'); % Nonlinear least squares exit flag
%
% Loop through Slices
%
for k = 1:nsl
%
   sl_time = zeros(nslt,1);            % Trigger time from files
   sll = int2str(k);                   % Slice number as letters
%
% Loop through Spin Lock Times
%
   for l = 1:nslt       % Loop through spin lock times
%
      n = nslt*(k-1)+l;
      fnam = fnams{id(n)};             % Filename for this spin lock time
      fnam = fullfile(dir_nam,fnam);
      fprintf(1,['\n Processing file:  ' strrep(fnam,'\','\\') ...
                 ', Slice:  ' sll ', Spin lock time:  ' ...
                 int2str(slt(l)) ' ms']);
%
% Load and Scale Image
%
      img = dicomread(fnam);
      info = dicominfo(fnam);
      sl_time(l) = info.TriggerTime;
      sl = single(info.RescaleSlope);
      offst = single(info.RescaleIntercept);     % Usually zero
      img = single(img);
      r = (img-offst)/sl;
      dat3d(:,:,l) = double(r);
   end
%
   fprintf(1,'\n');     % Line between slices
%
% Check Spin Lock Times
%
   if ~all(sl_time==slt)
     error(' *** ERROR in T1rho_map2:  Spin lock times not correct!');
   end
%
% Setup DICOM Header Information to Write Out to DICOM File
%
   info.SeriesDescription = ds;
   info.LargestImagePixelValue = maxint;
   info.RescaleSlope = rslope;
   info.RescaleIntercept = 0;
%
   info_res = info;
   info_res.RescaleSlope = ressl;
%
% Don't Fit Pixels with Small Initial Values
%
   for m = 1:nslt
      zvec = dat3d(:,:,m)<10^-m;       % Small (zero) values for this spin lock time
      zmat(:,:,k) = zmat(:,:,k)|zvec;  % Small (zero) values for this slice
   end
   vmat = repmat(~zmat(:,:,k),1,1,nslt);    % Valid values for this slice
%
% Get Log for Linear Least Squares
%
   dat3dl = double(vmat);
   dat3dl(vmat) = log(double(dat3d(vmat)));
%
% Calculate Linear Least Squares Solution by Rows
%
   for l = 1:npx
      row = squeeze(dat3d(l,:,:))';    % Row data
      rowl = squeeze(dat3dl(l,:,:))';  % Log row data
      p = xdat\rowl;
      T1rholls(l,:,k) = single(-1./p(1,:));   % Time constant in ms
      lls_amp(l,:,k) = single(p(2,:)); % Amplitude
      yp = exp(repmat(p(2,:),nslt,1)+xdat(:,1)*p(1,:));
      d1 = row-yp;
      d1 = sqrt(sum(d1.*d1));
      d2 = sqrt(sum(row.*row));
      lls_nres(l,:,k) = single(d1./d2);     % Normalized residuals
   end
%
% Write Linear Least Squares T1rho and Normalized Residuals to DICOM
% Files
%
   lidx = isinf(T1rholls(:,:,k));
   idx1 = find(lidx);
   T1rholls(idx1+(k-1)*npx2) = 0;      % Set Inf to zero
   dnam = ['T1rho_lls_' fs '_' sll '.dcm'];
   dicomwrite(uint16(T1rholls(:,:,k)/rslope),dnam,info, ...
              'CreateMode','Copy');
%
   rnam = ['lls_nres_' fs '_' sll '.dcm'];
   dicomwrite(uint16(lls_nres(:,:,k)/ressl),rnam,info_res, ...
              'CreateMode','Copy');
%
% Nonlinear Least Squares Exponential Fit to Get T1rho Values
%
   if ~nskip
     for l = 1:npx      % Rows
        parfor m = 1:npx               % Columns
           if vmat(l,m,1)              % Only calculate if valid data
             ydat = squeeze(dat3d(l,m,:));
             if T1rholls(l,m,k)>0&&T1rholls(l,m,k)<100     % Use linear least squares as initial parameters
               rp0 = double([lls_amp(l,m,k) T1rholls(l,m,k)]);
             else
               rp0 = [max(ydat) T1rho0];
             end
%
             [rp,rnorm,~,eflag] = lsqcurvefit(fun,rp0,slt,ydat, ...
                                              [],[],opt);
%
             if isinf(rp(2))
               T1rhonls(l,m,k) = 0;    % Set Inf to zero
             else
               T1rhonls(l,m,k) = single(rp(2));
             end
             nls_amp(l,m,k) = single(rp(1));
             nls_nres(l,m,k) = single(sqrt(rnorm)./norm(ydat));
             exit_flag(l,m,k) = int8(eflag);
           end
        end
     end
   end
%
% Write Nonlinear Least Squares T1rho and Normalized Residuals to DICOM
% Files
%
   dnam = ['T1rho_nls_' fs '_' sll '.dcm'];
   dicomwrite(uint16(T1rhonls(:,:,k)/rslope),dnam,info, ...
              'CreateMode','Copy');
%
   rnam = ['nls_nres_' fs '_' sll '.dcm'];
   dicomwrite(uint16(nls_nres(:,:,k)/ressl),rnam,info_res, ...
              'CreateMode','Copy');
%
% Plot Selective Slices
%
   if any(plt_sl==k)
     img = T1rholls(:,:,k);
     img(img<0) = 0;
     img(img>200) = 0;
     img(isnan(img)) = 0;
%
     figure;
     imagesc(img,[0 100]);
     axis off;
     axis image;
     title(['LLS T1rho for Slice ' sll],'FontSize',16,'FontWeight', ...
            'bold');
     colormap gray;
     colorbar;
     if iprt
       orient landscape;
       if k==plt_sl(1)
         print('-dpsc2','-r600','-fillpage',pnam);
       else
         print('-dpsc2','-r600','-fillpage','-append',pnam);
       end
     end
%
     figure;
     imagesc(lls_nres(:,:,k),[0 1]);
     axis off;
     axis image;
     title(['LLS Normalized Residuals for Slice ' sll],'FontSize', ...
            16,'FontWeight','bold');
     colormap gray;
     colorbar;
     if iprt
       orient landscape;
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
     img = T1rhonls(:,:,k);
     img(img<0) = 0;
     img(img>200) = 0;
     img(isnan(img)) = 0;
%
     figure;
     imagesc(img,[0 100]);
     axis off;
     axis image;
     title(['NLS T1rho for Slice ' sll],'FontSize',16,'FontWeight', ...
            'bold');
     colormap gray;
     colorbar;
     if iprt
       orient landscape;
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
     figure;
     imagesc(nls_nres(:,:,k),[0 1]);
     axis off;
     axis image;
     title(['NLS Normalized Residuals for Slice ' sll],'FontSize', ...
            16,'FontWeight','bold');
     colormap gray;
     colorbar;
     if iprt
       orient landscape;
       print('-dpsc2','-r600','-fillpage','-append',pnam);
     end
%
  end
%
end
%
% Save MAT File
%
mnam = ['T1rho_' fs '.mat'];
%
save(mnam,'id','id1','idx','fnams','nslt','slt','npx','T1rholls', ...
          'lls_amp','lls_nres','T1rhonls','nls_amp','nls_nres', ...
          'exit_flag','series_desc');
%
return
