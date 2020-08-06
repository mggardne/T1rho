%#######################################################################
%
%                    * T1rho SLice VERIFY Program *
%
%          M-File which plots ROIs on T1rho slices before and after the
%      digitized slice to verify the use of the correct slice.
%
%     NOTES:  1.  Run after running main plotting program
%            T1rho_sl_plt.m or loading T1rho MAT file.
%
%     RESULTS:  1.  Visually looks like the correct slices were used to
%               generate masks.
%
%     04-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Series Description and ROI Names For Legend Entries
%
fs = series_desc{idx};
fs = fs(8:end);         % Short series description w/o '3DMAPPS'
fs = strtrim(fs);       % File name description
%
rnams = {rois.name}';   % Get names of ROIs
legds = strrep(rnams,'.csv','');       % Create legend entries
%
% Loop through ROI Slices and Plot ROIs
%
lt = ['b.-'; 'g.-'; 'r.-'; 'c.-'; 'm.-'; 'y.-' ];  % Line color and type
%
nrsl = size(rsl,1);     % Number of slices with ROIs
%
for k = 1:nrsl
%
% Plot T1rho Slice Images
%
   slk = rsl(k);        % Slice number in spin lock images
   sl = irsl(k);        % Slice number in T1rho images
%
   for n = -1:1         % Check one slice before and after original slice
      figure;
      orient landscape;
      imagesc(T1rhonls(:,:,sl+n),[0 100]);
      colormap gray;
      axis image;
      axis off;
      title([fs ' T1rho Slice ' int2str(sl+n)],'FontSize',16, ...
            'FontWeight','bold');
      hold on;
%
      lh = gobjects(nrfiles,1);        % Line graphic handles
      idl = false(nrfiles,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
      for l = 1:nrfiles
         idxr = rois(l).slice==slk;
         if any(idxr)
           dat = cell2mat(rois(l).roi(idxr).data);
           lh(l) = plot(dat(:,1),dat(:,2),lt(l,:));
           idl(l) = true;
         end
      end
%
% Add Legends and Print Slice Plots
%
      legend(lh(idl),legds(idl));
%
   end                  % End of n loop - 1 slice before through 1 slice after
end                     % End of k loop - ROI slices loop
%
return