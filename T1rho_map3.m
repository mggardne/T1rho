%#######################################################################
%
%                       * T1rho_map 3 Program *
%
%          M-File which reads knee MRI data with different spinlock
%     times and fits an exponential to the data using nonlinear least
%     squares.  The T1rho map and normalized residuals are output to
%     DICOM format files.
%
%     NOTES:  1.  dicom_lst.m must have produced a MAT file with a list
%             of MRI DICOM files in the current directory.
%
%             2.  Must have M-File exp_fun1.m and application file
%             SeriesTable.mlapp in the current directory.
%             or path.
%
%     18-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Print Parameter
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots
%
% Check Spin Lock Times?
%
% slt_chk = true;         % Check spin lock times
slt_chk = false;        % Don't check spin lock times (Philips' TriggerTime is not correct.)
%
% Nonlinear Least Squares Parameter
%
nskip = false;          % Don't skip nonlinear calculations
% nskip = true;           % Skip nonlinear calculations
%
% Get DICOM File Names
%
matnam = 'dicom_lst.mat';
if exist(matnam,'file')
  load(matnam,'idv','sn','id','t0','afiles','splckt');
else
  error(' *** ERROR in T1rho_map3:  File dicom_lst.mat not found!');
end
%
% Let User Choose DICOM Series
%
ids = id(idv);          % Only series with spin lock times
nss = sum(idv);         % Number of series with spin lock times
%
t1 = t0(idv,[1:3 5 7 8 16:18 21:22]);  % Display a subset of whole table
t2 = false(nss,1);      % Column for choosing a series
t2 = table(t2,'VariableNames',{'Use for fitting?'});
t1 = [t1 t2];
%
clear t0 t2;
%
idc = 0;
%
while sum(idc)~=1       % Only one series at a time
%
     ht = SeriesTable(t1);
     uiwait(ht.UIFigure);
%
     idc = ht.UITable.Data{:,12};      % Series for T1rho calculations
%
     ht.delete;              % Delete figure with table
%
end
%
% Number of Spin Lock Times and Pixels for this Series
%
slt = eval([ '[' splckt{idc} ']' ])';
nslt = size(slt,1);     % Number of spin lock times
%
npx = t1.ImageColumns(idc);            % Number of pixels in X
npy = t1.ImageRows(idc);               % Number of pixels in Y
nps = npx*npy;          % Number of pixels in slice
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on');
%
T1rho0 = 80;            % Initial value for T1rho
fun = @exp_fun1;        % Exponential function
%
% Get Index to Files and Get File Names
%
ids = id(idv);
ids = ids(idc);         % Index to series
%
fnams = afiles{ids};
nf = size(fnams,1);     % Number of files
nsl = nf/nslt;          % Number of slices
%
plt_sl = floor(nsl/4);
plt_sl = plt_sl:plt_sl:3*plt_sl;       % Slices to plot
%
%  Get File Name Description and DICOM Header Parameters
%
sn1 = sn(idv);
sn1 = sn1(idc);
fs = ['S' int2str(sn1)];
ds = ['T1rho' fs];      % New DICOM description with series number
%
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
      fnam = fnams{n};  % Filename for this spin lock time
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
   if slt_chk
     if ~all(sl_time==slt)
       error(' *** ERROR in T1rho_map3:  Spin lock times not correct!');
     end
   end
%
% Setup DICOM Header Information to Write Out to DICOM File
%
   info.SeriesDescription = ds;
   info.LargestImagePixelValue = maxint;
   info.RescaleSlope = rslope;
   info.RescaleIntercept = 0;
   info.BodyPartExamined = 'KNEE';     % Must be upper case
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
   T1rholls(idx1+(k-1)*nps) = 0;      % Set Inf to zero
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
save(mnam,'id','idc','ids','idv','fnams','fs','nslt','slt','npx', ...
          'T1rholls','lls_amp','lls_nres','T1rhonls','nls_amp', ...
          'nls_nres','exit_flag','sn','sn1','t1');
%
return