%#######################################################################
%
%                      * T1rho AIR MASK Program *
%
%          M-File which creates a mask for areas of air in the slices
%      with ROIs and plots the ROIs.  The air mask is saved to the
%      MAT file.
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
% Load T1rho MAT File Data
%
mnam = 'T1rho_voxel0.5mmslice2.0.mat';
load(mnam,'npx','idx','irsl','series_desc','T1rhonls');
%
% Air Mask
%
nrsl = size(irsl,1);
maska = false(npx*npx,nrsl);
%
% Define Areas of Air in First Slice (Manually from previous plots)
%
xy1_1 = [1 1; 1 60; 25 1];             % Upper left triangle
xy1_2 = [1 245; 1 512; 105 512];       % Lower left triangle
%
% Create Mask for First Slice
%
mask1 = tri_mask(xy1_1,npx);           % Upper left triangle
mask2 = tri_mask(xy1_2,npx);           % Lower left triangle
maska(:,1) = mask1|mask2;
%
% Define Areas of Air in Second Slice (Manually from previous plots)
%
xy2_1 = NaN(3,2);                      % Blank (none)
xy2_2 = [1 245; 1 512; 125 512];       % Lower left triangle
%
% Create Mask for Second Slice
%
maska(:,2) = tri_mask(xy2_2,npx);
%
% Define Areas of Air in Third Slice (Manually from previous plots)
%
xy3_1 = [1 1; 1 40; 15 1];             % Upper left triangle
xy3_2 = [1 180; 1 512; 165 512];       % Lower left triangle
%
% Create Mask for Third Slice
%
mask1 = tri_mask(xy3_1,npx);           % Upper left triangle
mask2 = tri_mask(xy3_2,npx);           % Lower left triangle
maska(:,3) = mask1|mask2;
%
% Define Areas of Air in Slice 4 (Manually from previous plots)
%
xy4_1 = [1 1; 1 512; 25 512; 25 1];    % Left rectangle
xy4_2 = [1 135; 1 512; 185 512];       % Lower left triangle
%
mask1 = false(npx*npx,1);
mask1(1:25*512) = true;                % First 25 columns (rectangle)
mask2 = tri_mask(xy4_2,npx);           % Lower left triangle
maska(:,4) = mask1|mask2;
%
% Get ROI Coordinates into Arrays for Plotting
%
xp = NaN(2,2,4);
yp = NaN(2,2,4);
%
id1 = [2 3];            % Index to hypotenuse points for triangle 1
id1s = int2str(id1);
id2 = [1 3];            % Index to hypotenuse points for triangle 2
id2s = int2str(id2);
%
id4 = [3 4];            % Index to left edge of rectangle
id4s = int2str(id4);
%
for k = 1:nrsl
   l = int2str(k);
   if k==4
     id1s = id4s;
   end
   xp(:,:,k) = eval(['[xy' l '_1([' id1s '],1) xy' l '_2([' id2s, ...
                     '],1)]']);
   yp(:,:,k) = eval(['[xy' l '_1([' id1s '],2) xy' l '_2([' id2s, ...
                     '],2)]']);
end
%
% Loop Through the Slices and Plot the Air ROIs
%
fs = series_desc{idx};
fs = fs(8:end);         % Short series description w/o '3DMAPPS'
fs = strtrim(fs);       % File name description
%
anam = [fs '_AirROIs.ps'];            % Air ROIs print file name
%
amap = gray(256);       % Color map for air area
amap(1,:) = [0 0.8 0];                 % Air is green
%
for k = 1:nrsl
%
   sl = irsl(k);
%
   figure;
   orient landscape;
   imga = T1rhonls(:,:,sl);
   imagesc(imga,[0 100]);
   colormap gray;
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
   hold on;
   plot(xp(:,:,k),yp(:,:,k),'g-');
%
   if k==1
     print('-dpsc2','-r600','-fillpage',anam);
   else
     print('-dpsc2','-r600','-fillpage','-append',anam);
   end
%
   figure;
   orient landscape;
   imgm = imga;
   imgm(maska(:,k)) = -0.5;
   imagesc(imgm,[-0.5 100]);
   colormap(amap);
   axis image;
   axis off;
   title([fs ' T1rho Slice ' int2str(sl)],'FontSize',16, ...
         'FontWeight','bold');
%
   print('-dpsc2','-r600','-fillpage','-append',anam);
%
end
%
% Save Mask to MAT File
%
save(mnam,'maska','-append');
%
return