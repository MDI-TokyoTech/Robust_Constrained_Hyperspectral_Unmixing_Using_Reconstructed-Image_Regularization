function [abundanciesCollection] = getAbundanciesCollectionGaussianFields(n,ne,height,width,type,theta1,theta2)
% [abundanciesCollection] = getAbundanciesCollectionGaussianFields(n,ne,height,width,type,theta1,theta2)
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
% Returns a collection of synthetic samples of abundancy images by Gaussian Fields
% method. This work is based on:
%   Computations With Gaussian Random Fields
%   Kozintsev, Boris
%   1999, PhD Thesis 1999-3
%   ISR (Institute for System Research), University of Maryland
%
% Input:
%   - n: number of samples to build (> 0, defaults are 10).
%   - ne: number of images in each sample (>= 2, default is 5),.
%   - height,width: spatial dimensions of the images (> 0, defaults are 256,256).
%   - type: type of gaussian field used:
%       1 - spheric (default)
%		2 - exponential
%		3 - rational
%		4 - matern
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
%   - abundanciesCollection: a 4-dimensional matrix (n x width x height x ne).

%% Check parameters
% deafults
switch nargin
    case 0
        n = 10; ne = 5; height = 256; width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 1
        ne = 5; height = 256; width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 2
        height = 256; width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 3
        width = 256; type = 1; theta1 = 100; theta2 = 0;
    case 4
        type = 1; theta1 = 100; theta2 = 0;
    case 5
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
    case 6
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
% Check n
if n <= 0
	error('n has to be greater than zero');
end
% More checks made in getAbundanciesSampleGaussianField;

%% Build samples
% Initialize collection
abundanciesCollection = zeros(n,height,width,ne);
% Build
for i=1:n
    fprintf('Synthesizing sample: %g.\n',i);
    abundanciesCollection(i,:,:,:) = getAbundanciesSampleGaussianFields(ne, height, width, type, theta1, theta2);
end
fprintf('Collection completed.\n');
