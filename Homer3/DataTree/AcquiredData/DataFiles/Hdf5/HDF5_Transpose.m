function val = HDF5_Transpose(val, options)

if ~exist('options','var')
    options = '';
end

% Matlab stores contiguous muti-dimensional arrays in column-major order.
% HDF5 stores them in row-major order. We want to transpose the data to agree
% with the file format's storage order.
if (~isrow(val) && ~iscolumn(val))              || ...      % Matrices
     ~isempty(findstr('multidim', options))     || ...      % Force multi-dimensional even if vector
     ~isempty(findstr('2D', options))           || ...      % Force 2D even if vector
     ~isempty(findstr('3D', options))                       % Force 3D even if vector
    
    val = permute(val, ndims(val):-1:1);
    
elseif iscolumn(val) && ischar(val)
    
    val = permute(val, ndims(val):-1:1);
    
end

