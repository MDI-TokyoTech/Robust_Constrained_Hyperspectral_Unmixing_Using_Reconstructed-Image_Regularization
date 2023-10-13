function [syntheticImage,abundanciesGT] = getHyperspectralSyntheticImageLegendre(endmembers)
% [syntheticImage,abundanciesGT] = getHyperspectralSyntheticImageLegendre(endmembers)
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
%    You should have received a copy of the GNU General Public License
%    along with HYDRA.  If not, see <http://www.gnu.org/licenses/>.
%
% Return an hyperspectral synthetic image by Legendre polynomies method.
%
% Input:
%   - endmembers: the set of endmembers to use. The endmembers set must be
%   an e rows x p columns matrix, being 'e' the number of endmembers and 'p' the
%   number of bands of each endmember.
% Output:
%   - syntheticImage: a 3-dimensional matrix (w x h x p) representing the built synthetic hyperspectral image.
%   - abundanciesGT: a 3-dimensional matrix (w x h x e) representing the
%   abundancies images used to build the synthetic image.

%% Check parameters
if nargin < 1
    error('You must enter the set of endmembers to be used.');
end
% Check dimensions
[e p] = size(endmembers);
% Build default abundancies
fprintf('Building default abundancies sample\n');
abundanciesGT = getAbundanciesSampleLegendre(e);
% Check dimensions
[w h n] = size(abundanciesGT);
if e ~= n
    error('The number of endmembers and abundancies images must be the same.');
end

%% Synthetising
fprintf('Building image\n');
tic;
% Empty image
syntheticImage = double(zeros(w,h,p));
% Build the image
for i=1:w
    for j=1:h
        for k=1:e
            syntheticImage(i,j,:) = squeeze(syntheticImage(i,j,:))' + abundanciesGT(i,j,k)*endmembers(k,:);
        end
    end
end

t = toc;
fprintf('Synthetic hyperspectral image completed in %g seconds.\n',t);