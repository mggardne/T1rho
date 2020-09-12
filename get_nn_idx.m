function [idnn,idv] = get_nn_idx(idele,npx,ichk)
%GET_NN_IDX  Given a linear index into a two-dimensional (2D) square
%          matrix, returns a linear index matrix to the nine (9)
%          nearest neighbors.  The matrix has nine rows and the number
%          of columns equal the number of input points.
%
%          IDNN = GET_NN_IDX(IDELE,NPX) given the linear index into a
%          two-dimensional (2D) square matrix, IDELE, with dimensions
%          of NPX by NPX, returns a linear index matrix to the nine (9)
%          nearest neighbors, IDNN.  IDNN has dimensions of nine by 
%          the number of elements in the input index.
%
%          Edge elements which do not have nine nearest neighbors has
%          NaNs for indices to missing nearest neighbors.
%          
%          [IDNN,IDV] = GET_NN_IDX(IDELE,NPX) returns a logical matrix
%          of the same size as IDNN and is true for valid nearest
%          neighbors and false for missing nearest neighbors.
%
%          IDNN = GET_NN_IDX(IDELE,NPX,ICHK) if ICHK is true, the
%          nearest neighbors only includes elements in the input index.
%          This is useful for finding nearest neighbors within a region
%          of interest (ROI) within the 2D square matrix.  By default,
%          ICHK is false.
%
%          NOTES:  1.  Edge elements at the top of the 2D square matrix
%                  columns do not wrap around to reference elements at
%                  the bottom and bottom column elements do not wrap
%                  around and reference elements at the top of the
%                  column.
%
%                  2.  Similar to Matlab command conndef(2,'max'), but
%                  gets linear index to surrounding pixels instead of a
%                  matrix of ones.
%
%                  3.  Means and standard deviations of the nearest
%                  neighbors can be calculated using stat_nan.m.  See
%                  stat_nan.m.
%
%                  4.  Future version may extend this function to 2D
%                  rectangular matrices.
%
%          10-Sep-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<2
  error(' *** ERROR in GET_NN_IDX:  Two inputs are required!');
end
%
if nargin<3
  ichk = false;         % Nearest neighbors are not confined to input index
end
%
if isempty(ichk)
  ichk = false;         % Nearest neighbors are not confined to input index
end
%
% Check Size of 2D Matrix and Get Number of Elements in 2D Matrix
%
ndim = size(npx(:),1);
if ndim>2||ndim<1
  error([' *** ERROR in GET_NN_IDX:  Incorrect number of 2D ', ...
         'matrix dimensions!']);
end
%
if ndim==2
  npx = npx(1);         % Assumes 2D matrix is square
end
%
npx2 = npx*npx;         % Maximum index into 2D square matrix
%
% Get Input Linear Index as a Row Vector
%
idele = idele(:)';
npts = size(idele,2);
%
% Get Offsets to All Nearest Neighbors
%
% There are nine (9) nearest neighbors.  Similar to Matlab command
% conndef(2,'max'), but gets index offset to surrounding pixels in a
% two-dimensional (2D) image matrix with number of rows "npx".
%
inc = -1:1;             % Nearest neighbors including input index
%
ids = [npx npx npx]'*inc;
ids = ids+repmat(inc',1,3);
ids = ids(:);
%
ids = repmat(ids,1,npts);
%
% Add Element Index to Get Index to Element and Neighbors
%
idnn = repmat(idele,9,1)+ids;
%
% Check for Invalid Indices
%
% Elements before first column (<1) or after last column (>npx2).
%
idb1 = idnn(:)<1;
idb2 = idnn(:)>npx2;
%
if any(idb1)
  idnn(idb1) = NaN;
end
%
if any(idb2)
  idnn(idb2) = NaN;
end
%
% Check Top and Bottom Rows to Prevent Wrap Around
%
toprow = 1:npx:npx2;
itop = ismember(idele,toprow);
if any(itop)
  idnn(1:3:9,itop) = NaN;
end
%
botrow = npx:npx:npx2;
ibot = ismember(idele,botrow);
if any(ibot)
  idnn(3:3:9,ibot) = NaN;
end
%
% Use Only Nearest Neighbors Within Input Point IDs?
%
if ichk
%
  idnc = idnn(:);
  idnc = ~ismember(idnc,idele);
  idnn(idnc) = NaN;
%
end
%
% Output Logical Matrix to Valid Indices?
%
if nargout>1
  idv = ~isnan(idnn);
end
%
return