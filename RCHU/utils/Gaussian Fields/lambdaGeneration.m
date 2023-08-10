function [lambda] = lambdaGeneration(ndata, type, theta1, theta2)
% [lambda] = lambdaGeneration(ndata, type, theta1, theta2)
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

%% Generation
[lambda c] = generation(ndata, ndata, type, theta1, theta2);
min(lambda);
if min(lambda) < 0,
   warning('Some negative lambda value. New trial.');
   lambda = lambdaGeneration(ndata*2, type, theta1, theta2);
   min(lambda);
end
