function [z] = renormalize(z)
% [z] = renormalize(z1, z2, nh)
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

%% Renormalization
[n,m,nh]=size(z);
for i=1:n
	for j=1:m
		[vmax,dominante]=max(z(i,j,:));
		restosuma=0;
		for b=1:nh
			if b~=dominante
				restosuma=restosuma+z(i,j,b);
			end
		end
		for b=1:nh
			if b~=dominante
				if restosuma>0
					z(i,j,b)=(z(i,j,b)/restosuma)*(1-vmax);
				else
					z(i,j,b)=0;
				end
				
			end
		end
	end
end
