function roi = rd_roi5(filenam,ipx)
%RD_ROI5  Reads an OSIRIX ROI CSV file.
%
%         ROI = RD_ROI5(FILENAM) reads the OSIRIX ROI CSV file, FILENAM,
%         and returns the structure, ROI, with the names of the ROI in
%         field "name" and the X, Y and Z data points for each slice is
%         in the columns of cell arrays in field "data". The X, Y and Z
%         data is in ordered triplets.  Each triplet represent a point
%         from a slice in the ROI.  Each ROI is returned in a row in
%         the structure.  
%
%         ROI = RD_ROI5(FILENAM,IPX) if IPX is true (nonzero) the X and
%         Y pixel data points for each slice is in the columns of cell
%         arrays in field "data" in the structure ROI.  Also returns the
%         field "imageno" with the image (or slice) number.
%
%         NOTES:  1.  The data was collected from OSIRIX using the
%                 polygon tool.
%
%                 2.  If the ROI name does not start with an alphabetic
%                 letter, the letter "a" is prepended to the ROI name.
%
%         01-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in RD_ROI5:  An input file name is required!');
end
%
if (nargin<2)
  ipx = false;
end
%
if isempty(ipx)
  ipx = false;
end
%
% Open File and Read First Two Lines
%
fid = fopen(filenam,'rt');
lin = fgetl(fid);       % Read first line of headers
hdrs = textscan(lin,'%s','Delimiter',',');
hdrs = hdrs{1};
idn = find(startsWith(hdrs,'NumOfPoints'));
lin = fgetl(fid);       % Read first line of data
%
% Keep Reading Data Until End of File
%
while lin~=-1
     idx = strfind(lin,',');
     idx = [idx length(lin)+1];
     if ipx
       img_num = str2double(lin(1:idx(1)-1))+1;  % Image number
     end
     rnam = lin(idx(7)+1:idx(8)-1);    % Name of ROI
%
% Strip Any Quotes from ROI Name
%
     idq = rnam==''''|rnam=='"';       % Both single and double quotes
     idq = ~idq;
     rnam = rnam(idq);
%
% Check for Hyphens, Percents and Spaces
%
     htrap = strfind(rnam,'-');        % Trap for hyphens
     if ~isempty(htrap)
       rnam(htrap) = '_';              % Replace hyphens w/ underscore
     end
     htrap = strfind(rnam,'%');        % Trap for percents
     if ~isempty(htrap)
       rnam(htrap) = 'p';              % Replace percents w/ p's
     end
     htrap = strfind(rnam,' ');        % Trap for spaces
     if ~isempty(htrap)
       rnam(htrap) = '_';              % Replace spaces w/ underscore
     end
%
% Check that First Character is a Letter
%
     if ~isletter(rnam(1))
       rnam = ['a' rnam];              % Prepend an "a" if not a character
     end
%
% Number of Data Points for this ROI
%
     npts = eval(lin(idx(idn-1)+1:idx(idn)-1));
     idp = idn+(0:5:(npts-1)*5)';
     if ipx
       mat = zeros(npts,2);
     else
       mat = zeros(npts,3);
     end
%
% Get Point Data
%
     for k = 1:npts
        if ipx
          mat(k,:) = eval(['[' lin(idx(idp(k)+3)+1:idx(idp(k)+5)-1) ']']);
        else
          mat(k,:) = eval(['[' lin(idx(idp(k))+1:idx(idp(k)+3)-1) ']']);
        end
     end
%
% Save Data for Each ROI
%
     if exist(rnam,'var')
       eval(['nslice = size(' rnam ',2);']);
       eval([rnam '{' int2str(nslice+1) '} = mat;']);
       if ipx
         eval([rnam '{2,' int2str(nslice+1) '} = img_num;']);
       end
     else
       eval([rnam '{1} = mat;']);
       if ipx
         eval([rnam '{2,1} = img_num;']);
       end
     end
     lin = fgetl(fid);
end
%
% Close File
%
fclose(fid);
%
% Clear Workspace
%
clear ans fid filenam hdrs htrap idn idq lin k idx img_num ipx nslice ...
      rnam npts idp mat;
%
% Get ROI
%
rnam = who;             % Names of ROI
nvar = size(rnam,1);    % Number of ROI
ipx = eval(['size(' rnam{1} ',1)'])>1;
%
% Get Data into a Cell
%
data = cell(nvar,1);
if ipx
  img_num = cell(nvar,1);
end
%
for k = 1:nvar
  data{k} = eval([ rnam{k} '(1,:)']);
  if ipx
    img_num{k} = eval([ rnam{k} '{2,:}']);
  end
end
%
% Put Information into a Structure
%
roi = struct('name',rnam,'data',data);
if ipx
  [roi.imageno] = deal(img_num{:});
end
%
return