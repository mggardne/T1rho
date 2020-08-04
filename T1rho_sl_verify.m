%
% Plots ROIs on T1rho Slices Before and After Digitized Slice to
% Verify Correct Slice
%
% Note:  Run after running main plotting program T1rho_sl_plt.m or
% loading T1rho MAT file.
%
% Result:  Visually looks like correct slices were used to generate
% masks.
%
% 04-Aug-2020 * Mack Gardner-Morse
%

%
% Loop through ROI Slices and Plot ROIs
%
for k = 1:nrsl
%
% Plot T1rho Slice Images
%
   slk = rsl(k);        % Slice number in spin lock images
   sl = irsl(k);        % Slice number in T1rho images
%
   for m = -1:1         % Check one slice before and after original slice
   figure;
   orient landscape;
   imagesc(T1rhonls(:,:,sl+m),[0 100]);
   colormap gray;
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl+m)],'FontSize',16, ...
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
   legend(lh(idl),legds(idl));
   end                  % End of m loop - 1 slice before through 1 slice after
end
