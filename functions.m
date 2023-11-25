% WRITTEN BY AHMET KEMAL CABBAR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean (Average) Filtering with the user-entered kernel size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputImage = meanFilter(inputImage, kernelSize)
    [rows, cols] = size(inputImage);
    outputImage = zeros(rows, cols);

    % Ensure the kernel size is odd
    if mod(kernelSize, 2) == 0
        error('Kernel size must be odd.');
    end

    halfKernel = (kernelSize - 1) / 2;

    for i = 1:rows
        for j = 1:cols
            % Get the neighborhood
            neighborhood = inputImage(max(1, i - halfKernel):min(rows, i + halfKernel), ...
                                      max(1, j - halfKernel):min(cols, j + halfKernel));

            % Compute the mean of the neighborhood
            outputImage(i, j) = mean(neighborhood(:));
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Median Filtering with the user-entered kernel size,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputImage = medianFilter(inputImage, filterSize)
    
    % Get the dimensions of the input image
    [rows, cols] = size(inputImage);
    
    % Initialize the output image with zeros
    outputImage = zeros(rows, cols);
    
    % Calculate the half-size of the filter
    halfSize = floor(filterSize / 2);
    
    for i = 1:rows
        for j = 1:cols
            % Extract the neighborhood of the current pixel
            neighborhood = inputImage(max(i - halfSize, 1):min(i + halfSize, rows), ...
                                      max(j - halfSize, 1):min(j + halfSize, cols));
            
            % Apply median filtering to the neighborhood
            outputImage(i, j) = median(neighborhood(:));
        end
    end

    % it didn't work when i didnot the data type convertion
    outputImage = uint8(outputImage);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Median Filtering with the user-entered maximum window size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputImage = adaptiveMedianFilter(inputImage, maxWindowSize)
    [rows, cols] = size(inputImage);
    outputImage = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols
            % Get the current pixel and its neighborhood
            pixelValue = inputImage(i, j);
            windowSize = 3; % Initial window size

            while windowSize <= maxWindowSize
                % Get the neighborhood
                neighborhood = inputImage(max(1, i - (windowSize-1)/2):min(rows, i + (windowSize-1)/2), ...
                                          max(1, j - (windowSize-1)/2):min(cols, j + (windowSize-1)/2));

                % Check if the pixel value is the median within the neighborhood
                if pixelValue == median(neighborhood(:))
                    outputImage(i, j) = pixelValue;
                    break;
                else
                    % Increase the window size and try again
                    windowSize = windowSize + 2;
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Filtering (Smoothing) with the given standard derivation value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sonuc = GaussianKernel(std)

    size_number= std*6;
    
    %check the number
    %if it is even, add1
    if mod(size_number,2) ==0
        size_number = size_number + 1;
    end
    
    [X, Y] = meshgrid(-(size_number-1)/2:(size_number-1)/2, -(size_number-1)/2:(size_number-1)/2);
    
    firstterm=1/(2*pi*std^2);
    secondterm=exp(-1*(X.^2+Y.^2)/(2*std^2));
    
    sonuc=firstterm*secondterm;
    
    % Normalize the kernel
    sonuc = sonuc / sum(sonuc(:));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Additive Uniform Noise with user parameters,
% - Additive Gaussian Noise with user parameters,
% - Additive Salt&Pepper Noise with user parameters,
% - Additive LogNormal Noise with user parameters,
% - Additive Rayleigh Noise with user parameters,
% - Additive Exponential Noise with user parameters,
% - Additive Erlang Noise with user parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = imnoise2(type, M, N, a, b)
%IMNOISE2 Generates an array of random numbers with specified PDF.
%   R = IMNOISE2(TYPE, M, N, A, B) generates an array, R, of size
%   M-by-N, whose elements are random numbers of the specified TYPE
%   with parameters A and B. If only TYPE is included in the
%   input argument list, a  single random number of the specified
%   TYPE and default parameters shown below is generated. If only
%   TYPE, M, and N are provided, the default parameters shown below
%   are used.  If M = N = 1, IMNOISE2 generates a single random
%   number of the specified TYPE and parameters A and B.
%
%   Valid values for TYPE and parameters A and B are:
% 
%   'uniform'       Uniform random numbers in the interval (A, B).
%                   The default values are (0, 1).  
%   'gaussian'      Gaussian random numbers with mean A and standard
%                   deviation B. The default values are A = 0, B = 1.
%   'salt & pepper' Salt and pepper numbers of amplitude 0 with
%                   probability Pa = A, and amplitude 1 with
%                   probability Pb = B. The default values are Pa =
%                   Pb = A = B = 0.05.  Note that the noise has
%                   values 0 (with probability Pa = A) and 1 (with
%                   probability Pb = B), so scaling is necessary if
%                   values other than 0 and 1 are required. The noise
%                   matrix R is assigned three values. If R(x, y) =
%                   0, the noise at (x, y) is pepper (black).  If
%                   R(x, y) = 1, the noise at (x, y) is salt
%                   (white). If R(x, y) = 0.5, there is no noise
%                   assigned to coordinates (x, y). 
%   'lognormal'     Lognormal numbers with offset A and shape
%                   parameter B. The defaults are A = 1 and B =
%                   0.25.
%   'rayleigh'      Rayleigh noise with parameters A and B. The
%                   default values are A = 0 and B = 1. 
%   'exponential'   Exponential random numbers with parameter A.  The
%                   default is A = 1.
%   'erlang'        Erlang (gamma) random numbers with parameters A
%                   and B.  B must be a positive integer. The
%                   defaults are A = 2 and B = 5.  Erlang random
%                   numbers are approximated as the sum of B
%                   exponential random numbers.

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.5 $  $Date: 2003/10/12 23:37:29 $

% Set default values.
if nargin == 1
   a = 0; b = 1;
   M = 1; N = 1;
elseif nargin == 3
   a = 0; b = 1;
end

% Begin processing. Use lower(type) to protect against input being
% capitalized. 
switch lower(type)
case 'uniform'
   R = a + (b - a)*rand(M, N);
case 'gaussian'
   R = a + b*randn(M, N);
case 'salt & pepper'
   if nargin <= 3
      a = 0.05; b = 0.05;
   end
   % Check to make sure that Pa + Pb is not > 1.
   if (a + b) > 1
      error('The sum Pa + Pb must not exceed 1.')
   end
   R(1:M, 1:N) = 0.5;
   % Generate an M-by-N array of uniformly-distributed random numbers
   % in the range (0, 1). Then, Pa*(M*N) of them will have values <=
   % a. The coordinates of these points we call 0 (pepper
   % noise). Similarly, Pb*(M*N) points will have values in the range
   % > a & <= (a + b).  These we call 1 (salt noise). 
   X = rand(M, N);
   c = find(X <= a);
   R(c) = 0;
   u = a + b;
   c = find(X > a & X <= u);
   R(c) = 1;
case 'lognormal'
   if nargin <= 3
      a = 1; b = 0.25;
   end
   R = a*exp(b*randn(M, N));
case 'rayleigh'
   R = a + (-b*log(1 - rand(M, N))).^0.5;
case 'exponential'
   if nargin <= 3
      a = 1;
   end
   if a <= 0
      error('Parameter a must be positive for exponential type.')
   end
   k = -1/a;
   R = k*log(1 - rand(M, N));
case 'erlang'
   if nargin <= 3
      a = 2; b = 5;
   end
   if (b ~= round(b) | b <= 0)
      error('Param b must be a positive integer for Erlang.')
   end
   k = -1/a;
   R = zeros(M, N);
   for j = 1:b         
      R = R + k*log(1 - rand(M, N));
   end
otherwise
   error('Unknown distribution type.')
end
end