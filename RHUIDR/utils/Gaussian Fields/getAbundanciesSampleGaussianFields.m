function [abundanciesSample] = getAbundanciesSampleGaussianFields(ne, height, width, type, theta1, theta2)
% [abundanciesSample] = getAbundanciesSampleGaussianFields(ne, height, width, type, theta1, theta2)
%
%    Copyright 2007,2008 Grupo Inteligencia Computacional, Universidad del País Vasco / Euskal Herriko Unibertsitatea (UPV/EHU)
%
%    Website: http://www.ehu.es/ccwintco
%
%    This file is part of the HYperspectral Data Retrieval and Analysis tools (HYDRA) program.
%
%    HYDRA is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    HYDRA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with HYDRA.  If not, see <http://www.gnu.org/licenses/>.
%
% Returns a synthetic sample of abundancy images by Gaussian Fields
% method. This work is based on:
%   Computations With Gaussian Random Fields
%   Kozintsev, Boris
%   1999, PhD Thesis 1999-3
%   ISR (Institute for System Research), University of Maryland
%
% Input:
%   - ne: number of images in each sample (>= 2, default is 5),.
%   - height,width: spatial dimensions of the images (> 0, defaults are 256,256).
%   - type: type of gaussian field used:
%       1 - spheric (default)
%		2 - exponential
%		3 - rational
%		4 - matern
%           1,4両方試してみる
%   - theta1,theta2: parameters for the gaussian field type
%       - spheric (defaults):    
%                   theta1 > 0 (usually between 1 and 200) 
%			        theta2 = 0 (don't used)
%		- exponential:  
%                   theta1 = (0,1)
%                   theta2 = (0, 2]
%				Usually:
%					theta1 = 0.5, 0.7, 0.9
%					theta2 = 0.8, 1.3, 1.9
%
%		- rational:     
%                   theta1, theta2 > 0
%				Usually:
%					theta1 = 1.2, 1.6, 1.85
%					theta2 = 0.85, 1.5, 2 	
%		- matern: 	
%                   theta1, theta2 > 0
%				Usually:
%					theta1 = 5, 10, 20
%					theta2 = 0.5, 1, 2
%
% Output:
%   - abundanciesSample: a 3-dimensional matrix (width x height x ne).

%% Check parameters
% defaults
switch nargin
    case 0
        ne = 5; height = 256; width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 1
        height = 256; width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 2
        width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 3
        type = 1; theta1 = 100; theta2 = 0;
    case 4
        switch type
            case 1
                theta1= 100; theta2 = 0;
            case 2
                theta1= 0.7; theta2 = 1.3;
            case 3
                theta1= 1.6; theta2 = 1.5;
            case 4
                theta1= 10; theta2 = 1;
        end
    case 5
        switch type
            case 1
                theta2 = 0;
            case 2
                theta2 = 1.3;
            case 3
                theta2 = 1.5;
            case 4
                theta2 = 1;
        end
end
% check
if height <= 0 || width <= 0
	error('width and height has to be greater than zero\n');
end
if ne < 2
	error('ne has to be equal or greater than 2');
end
if type < 1 || type > 4
	error('gaussian field type is incorrect');
end
switch type
	case 1
        if (theta1 < 0)
        	warning('Model parameters are incorrect, defaults used.')
            theta1 = 100; theta2 = 0;
        end
	case 2
        if theta1 <= 0 || theta1 >= 1 || theta2 <= 0 || theta2 > 2
        	warning('Model parameters are incorrect, defaults used.')
        	theta1= 0.7; theta2 = 1.3;
        end
	case 3
        if theta1 <= 0 || theta2 <= 0
        	warning('Model parameters are incorrect, defaults used.')
        	theta1= 1.6; theta2 = 1.5;
        end
	case 4
        if theta1 <= 0 || theta2 <= 0
        	warning('Model parameters are incorrect, defaults used.')
        	theta1= 10; theta2 = 1;
        end
end

%% Build the sample
tic;
% Get the fewer power of 2 greater than width and height
nn = 2^nextpow2(max(height,width));
% Generate lambda
lambda = lambdaGeneration(nn, type, theta1, theta2);
% Generate abundancies images
[z1 z2] = getSamples(lambda, height,width, ne);
% Normalize abundancies images
[abundanciesSample zn2] = normalization(z1, z2, ne);
t = toc;
fprintf('Sample completed in %g seconds.\n',t);
