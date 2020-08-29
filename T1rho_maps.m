%#######################################################################
%
%                   * T1rho Map Sensitivity Program *
%
%          M-File which reads knee MRI data with different spinlock
%     times and fits an exponential to the cartilage data using
%     nonlinear least squares.  The fits are only done in areas with
%     cartilage on slices with cartilage regions of interest.  The fits
%     are done with different combinations of spin lock times to assess
%     the sensitivity of the results to the number of spin lock times.
%     Twelve random points on each slice are plotted with the different
%     fits.  The T1rho maps and square root of the sum of squared
%     errors are saved to Matlab MAT files.  Slice and overall
%     statistics are written to a MS-Excel spreadsheet.
%
%     NOTES:  1.  Must have M-files exp_fun1.m and ncombi.m in the
%             current directory or path.
%
%     26-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Print Parameter
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots
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
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on');
%
T1rho0 = 80;            % Initial value for T1rho
fun = @exp_fun1;        % Exponential function
%
%  Line Color/Type and Set Plot File Name
%
lnctyp = ['r-'; 'g-'; 'b-'];         % Line color and type
%
if iprt
  pnam = ['T1rho_maps_' fs '.ps'];     % Plot file name
end
%
% Combine Masks to Get Mask for All Cartilage on a Slice
%
mask1 = squeeze(maskf(:,1,:)|maskp(:,1,:)|maskt(:,1,:));   % Layer 1
mask2 = squeeze(maskf(:,2,:)|maskp(:,2,:)|maskt(:,2,:));   % Layer 2
maskc = mask1|mask2;    % All cartilage mask
%
% Get Combinations of Spin Lock Times
% Always use the two end spin lock times at 0 ms and 80 ms
% (index 1 and nslt)
%
nn = nslt-1;            % Number of spin lock times starting with two
ncomb = cell(nn,1);
nsz = zeros(nn,1);
v = 2:nn;               % Index to middle spin lock times
%
for k = 0:nn-1
   a = nchoosek(v,k);
   nsz(k+1) = size(a,1);
   ncomb{k+1} = [ones(nsz(k+1),1) a repmat(nslt,nsz(k+1),1)];
end
%
nt = sum(nsz);          % Number of trials
nszc = [0; cumsum(nsz(1:nn-1))];       % Index for trial counter
%
% Put the Different Spin Lock Time Combinations into Columns for Output
%
trials = NaN(nt,nslt);
%
for l = 1:nn
   for m = 1:nsz(l)
      idr = nszc(l)+m;            % Index to rows
      trials(idr,1:l+1) = slt(ncomb{l}(m,:))';
   end
end
%
% Set Up Arrays for Loop
%
xdat = [slt ones(nslt,1)];             % Spin lock times in ms and exponential amplitude
%
np = zeros(nrsl,1);     % Number of cartilage pixels per slice
%
T1rhonlss = cell(nrsl,1);              % Nonlinear least squares time constant
nlss_amp = cell(nrsl,1);               % Nonlinear least squares amplitude
nlss_sse = cell(nrsl,1);               % Nonlinear least squares square root of sum of squared errors
exit_flags = cell(nrsl,1);             % Nonlinear least squares exit flag
%
stats = cell(1,nrsl);   % Curvefit statistics
%
xnam = ['T1rho_maps_' fs '.xlsx'];     % MS-Excel File
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
% Get Log for Linear Least Squares
%
   dat2dl = log(dat2d);
%
% Calculate Linear Least Squares Solution
% Least Squares Parameters - r
%
   r = xdat\dat2dl;     % Exponent in first row and amplitude in second row
   r(1,:) = -1./r(1,:); % Time constant in ms
   r = flipud(r);       % Amplitude in first row and T1rho in second row
%
% Initialize Slice Arrays
%
   t1r = zeros(npk,nt,'single');       % T1rho (time constant)
   amp = zeros(npk,nt,'single');       % Amplitude
   sse = zeros(npk,nt,'single');       % Square root of sum of squared errors
   flags = zeros(npk,nt,'uint8');      % Sum of squared errors
%
% Nonlinear Least Squares Exponential Fit to Get T1rho Values
%
   for l = 1:nn         % Number of spin lock times to use in fit minus one
      for m = 1:nsz(l)  % Number of combinations with this number of spin lock times
%
         idx = nszc(l)+m;              % Index to columns
         nc = ncomb{l}(m,:);           % Index to spin lock time combinations
         slto = slt(nc);               % Reduced number of spin lock times
%
         parfor o = 1:npk              % Number of pixels
%
           ro = r(:,o);
           dato = dat2d(:,o);
           ydat = dato(nc);            % Only spin lock times in "nc"
           if ro(2)>0&&ro(2)<100       % Use linear least squares as initial parameters
             rp0 = ro;
           else
             rp0 = [max(ydat); T1rho0];
           end
%
           [rp,~,~,eflag] = lsqcurvefit(fun,rp0,slto,ydat,[],[],opt);
%
           if isinf(rp(2))
             t1r(o,idx) = 0;           % Set Inf to zero
           else
             t1r(o,idx) = single(rp(2));
           end
%
           amp(o,idx) = single(rp(1));
           d = exp_fun1(rp,slt)-dato;  % Use all spin lock times
           sse(o,idx) = single(sqrt(d'*d));
           flags(o,idx) = int8(eflag);
%
         end
      end
   end
%
% Put Nonlinear Least Squares T1rho and Errors into Cell Arrays
%
   T1rhonlss{k} = t1r;
   nlss_amp{k} = amp;
   nlss_sse{k} = sse;
   exit_flags{k} = flags;
%
% Get and Write Statistics to MS-Excel Spreadsheet
%
   npt = repmat(np(k),nt,1);
   np0 = sum(t1r>0&t1r<100)';
   np20 = sum(t1r>20&t1r<100)';
   t1rm = mean(t1r)';
   t1rmd = median(t1r)';
   t1rsd = std(t1r)';
   ssem = mean(sse)';
   ssemd = median(sse)';
   ssesd = std(sse)';
   stats{k} = [npt,np0,np20,t1rm,t1rmd,t1rsd,ssem,ssemd,ssesd];
%
% Plot Select Pixels
%
   for l = 1:12
      o = ipp(l);
      if l==1||l==7
       hf = figure;
       orient tall;
      end
      if l>6
        subplot(3,2,l-6);
      else
        subplot(3,2,l);
      end
%
      y = dat2d(:,o);
      plot(slt,y,'bo:','LineWidth',1); % Raw data for all spin lock times
      hold on;
      ic = cell(3,1);
      for m = 3:5       % Use from three to five spin lock times
         i1 = m-2;
         idcs = nszc(m-1)+1:nszc(m);
         [~,idmn] = min(sse(o,idcs)); % Get best fit
         idmn = idcs(idmn);
%
         [l1,m1] = ncombi(idmn,nszc);
         ic{i1} = int2str(slt(ncomb{l1}(m1,:))');
         ic{i1} = [ic{i1} ',rsse=' sprintf('%.3f',sse(o,idmn))];
%
         yf = exp_fun1([amp(o,idmn); t1r(o,idmn)],slt);
         plot(slt,yf,lnctyp(i1,:),'LineWidth',1);     % Curvefits
      end
%
      legend(['Raw Data'; ic],'FontSize',8,'Location','northeast');
      if any(l==[5 6 11 12])
        xlabel('Spin Lock Time (ms)','FontSize',10,'FontWeight', ...
               'bold');
      end
      if any(l==1:2:11)
        ylabel('T1 Signal','FontSize',10,'FontWeight','bold');
      end
      title(['Pixel ' ipps(l,:)],'FontSize',12,'FontWeight','bold');
      if l==6||l==12
        sgtitle([fs ' Slice ' int2str(irsl(k))],'FontSize',15, ...
                 'FontWeight','bold');
      end
      if iprt
        if k==1&&l==6
          print('-dpsc2','-r600','-fillpage',pnam);
        elseif l==6||l==12
          print('-dpsc2','-r600','-fillpage','-append',pnam);
        end
      end
%
   end
%
end
%
% Get Overall Statistics
%
t1r = cell2mat(T1rhonlss);             % Combine all slices
sse = cell2mat(nlss_sse);              % Combine all slices
%
npt = repmat(sum(np),nt,1);
np0 = sum(t1r>0&t1r<100)';
np20 = sum(t1r>20&t1r<100)';
t1rm = mean(t1r)';
t1rmd = median(t1r)';
t1rsd = std(t1r)';
ssem = mean(sse)';
ssemd = median(sse)';
ssesd = std(sse)';
%
% Create Table and Write Mean Values to MS-Excel Spreadsheet
%
spltc = cellstr([repmat('Spin Lock Time ',nslt,1) int2str((1:nslt)')])';
stato = {'# of Pixels Slice X','# of Valid Pixels 0-100 Slice X', ...
         '# of Valid Pixels 20-100 Slice X','Mean T1rho Slice X', ...
         'Median T1rho Slice X','SD T1rho Slice X', ...
         'Mean RSSE Slice X','Median RSSE Slice X','SD RSSE Slice X'};
%
nlab = size(stato,2);   % Number of labels per slice
nlt = nlab*nrsl;        % Total number of labels for all slices
statc = cell(1,nlt);
for k = 1:nrsl
   idlx = k*nlab;
   idlx = idlx-nlab+1:idlx;
   statc(idlx) = strrep(stato,'X',int2str(irsl(k)));
end
%
stato = strrep(stato,' Slice X','');
%
varnams = [spltc statc stato];
%
ts = table([trials,cell2mat(stats),npt,np0,np20,t1rm,t1rmd,t1rsd, ...
           ssem,ssemd,ssesd]);
ts = splitvars(ts);
ts.Properties.VariableNames = varnams;
%
writetable(ts,xnam);
%
% Plot of Mean Statistics by Number of Spin Lock Times
%
t1rr = zeros(nn,2);
sser = zeros(nn,2);
%
for l = 1:6
   if l==nn
     idrr = nszc(l)+1:nt;
   else
     idrr = nszc(l)+1:nszc(l+1);
   end
%
   t1rr(l,1) = mean(t1rm(idrr));
   t1rr(l,2) = mean(t1rmd(idrr));
   sser(l,1) = mean(ssem(idrr));
   sser(l,2) = mean(ssemd(idrr));
end
%
figure;
orient landscape;
%
x = (2:nslt)';
[ha,h1,h2] = plotyy(x,t1rr,x,sser);
%
h = [h1; h2];           % Handles for all four lines
set(h,'LineWidth',2);
clrs = get(h,'Color');
%
for k = 1:4
   set(h(k),'MarkerFaceColor',clrs{k});
end
%
set(ha,'FontSize',12,'FontWeight','bold','XTick',2:nslt);
xlabel('Number of Spin Lock Times','FontSize',16,'FontWeight','bold');
ylabel(ha(2),'RSSE (ms)','FontSize',16,'FontWeight','bold');
ylabel('T1\rho (ms)','FontSize',16,'FontWeight','bold');
grid on;
set(h1,'Marker','o','MarkerSize',7);
set(h2,'Marker','s','MarkerSize',8);
title('Mean Fit Statistics','FontSize',24,'FontWeight','bold');
legend(h,{'Mean T1rho';'Median T1rho';'Mean RSSE';'Median RSSE'}, ...
          'FontSize',16,'Location','northeast');
%
if iprt
  pnam2 = ['T1rho_maps2_' fs '.ps'];
  print('-dpsc2','-r600','-fillpage',pnam2);
end
%
% Save MAT File
%
smnam = ['T1rho_maps_' fs '.mat'];
%
save(smnam,'mask1','mask2','maskc','ncomb','nn','np','nrsl','nsz', ...
           'nszc','nt','T1rhonlss','nlss_amp','nlss_sse', ...
           'exit_flags','statc','stato','stats','trials','ts');
%
return