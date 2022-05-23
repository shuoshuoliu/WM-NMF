function mis = im2mis(varargin)
%IM2MIS Convert image to Java MemoryImageSource.
%
%   MIS = IM2MIS(I) converts the intensity image I to a Java
%   MemoryImageSource.
%
%   MIS = IM2MIS(X,MAP) converts the indexed image X with colormap MAP to a
%   Java MemoryImageSource.
%
%   MIS = IM2MIS(RGB) converts the RGB image RGB to a Java
%   MemoryImageSource.

%   Copyright 2001-2004 The MathWorks, Inc.  


%   This is really equivalent to images\images\im2mis.m Rev 1.3. 

%   Input-output specs
%   ------------------ 
%   I:    2-D, real, full matrix
%         uint8, uint16, or double
%         logical ok but ignored
%
%   RGB:  3-D, real, full matrix
%         size(RGB,3)==3
%         uint8, uint16, or double
%         logical ok but ignored
%
%   X:    2-D, real, full matrix
%         uint8 or double
%         if isa(X,'uint8'): X <= size(MAP,1)-1
%         if isa(X,'double'): 1 <= X <= size(MAP,1)
%         logical ok but ignored
%         
%   MAP:  2-D, real, full matrix
%         size(MAP,1) <= 256
%         size(MAP,2) == 3
%         double
%         logical ok but ignored
%
%   MIS:  Java MemoryImageSource

[img,map,method] = ParseInputs(varargin{:});

% Assign function according to method
switch method
  case 'intensity'
    mis = im2mis_intensity(img);
  case 'rgb'
    mis = im2mis_rgb(img);
  case 'indexed'
    mis = im2mis_indexed(img,map);
end    


%----------------------------------------------------
function mis = im2mis_intensity(I)

mis = im2mis_indexed(I,gray(256));


%----------------------------------------------------
function mis = im2mis_rgb(RGB)

mis = im2mis_packed(RGB(:,:,1),RGB(:,:,2),RGB(:,:,3));


%----------------------------------------------------
function mis = im2mis_packed(red,green,blue)

[mrows,ncols,~] = size(red);
alpha = 255*ones(mrows,ncols);
packed = bitshift(uint32(alpha),24);
packed = bitor(packed,bitshift(uint32(red),16));
packed = bitor(packed,bitshift(uint32(green),8));
packed = bitor(packed,uint32(blue));
pixels = packed';
mis = java.awt.image.MemoryImageSource(ncols,mrows,pixels(:),0,ncols);


%----------------------------------------------------
function mis = im2mis_indexed(x,map)

[mrows,ncols] = size(x);
map8 = uint8(round(map*255)); % convert color map to uint8
% Instantiate a ColorModel with 8 bits of depth
cm = java.awt.image.IndexColorModel(8,size(map8,1),map8(:,1),map8(:,2),map8(:,3));
xt = x';
mis = java.awt.image.MemoryImageSource(ncols,mrows,cm,xt(:),0,ncols);


%-------------------------------
% Function  ParseInputs
%
function [img, map, method] = ParseInputs(varargin)

% defaults
img = [];
map = [];
method = 'intensity';

if nargin<1
    return
end

img = varargin{1};

if ~isnumeric(img) || ~isreal(img) || issparse(img) 
    error(message('stats:im2mis:RealNotSparse'));
end


switch nargin
  case 1
    % figure out if intensity or RGB
    if ndims(img) == 2
        method = 'intensity';
    elseif ndims(img)==3 && size(img,3)==3
        method = 'rgb';
    else
        error(message('stats:im2mis:IntensityOrOther'));
    end
    
    % Convert to uint8.
    if isa(img,'double')
        img = uint8(img * 255 + 0.5);
        
    elseif isa(img,'uint16')
        img = uint8(bitshift(img, -8));
        
    elseif isa(img, 'uint8')
        % Nothing to do.
        
    else
        error(message('stats:im2mis:BadImage'));
    end
    
  case 2
    
    % indexed image
    method = 'indexed';
    map = varargin{2};
    
    % validate map
    if ~isnumeric(map) || ~isreal(map) || issparse(map) || ~isa(map,'double') 
        error(message('stats:im2mis:MapRealDouble'));
    end

    if size(map,2) ~= 3
        error(message('stats:im2mis:MapColumns'));
    end
    
    ncolors = size(map,1);
    if ncolors > 256
        error(message('stats:im2mis:MapColors'));
    end
    
    % validate img 
    if ndims(img) ~= 2
        error(message('stats:im2mis:MapTwoDim'));
    end
    
    if isa(img,'uint8')
        if max(img(:)) > ncolors-1
            error(message('stats:im2mis:MapLargeIndex'));
        end            
    elseif isa(img,'double')
        if max(img(:)) > ncolors
            error(message('stats:im2mis:MapLargeIndex'));
        end            
        if min(img(:)) < 1
            error(message('stats:im2mis:MapSmallIndex'));
        end
        
        img = uint8(img - 1);
    else
        error(message('stats:im2mis:IntOrDouble'));
    end
    
  otherwise
    error(message('stats:im2mis:Internal'));
end
