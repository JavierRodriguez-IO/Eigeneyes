function [dataZ_eigeneyes_pre] = eigeneyes_image_registration (origin, dataZ_eigeneyes_pre, Num_eyes)
% This script applies the registration of our data. 
% It places the (0,0) on the retina, on the cornea or on the anterior lens (by default).
% Parameters:
%   - origin: Indicates when you want to put the (0,0) in your coordinate
%             system. It can be 'retina', 'apex_cornea', or ' '.
%   - dataZ_eigeneyes_pre: Original Z-axis coordinates
%   - Num_eyes: Number of eyes in the database.
% Results:
%   - dataZ_eigeneyes_pre: New dataZ with the registration ((0,0)) in the first
%                          parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eigeneyes: A statistical model of whole-eye geometry and its applications
% Copyright (C) 2026 Estela Vadillo, Álvaro de la Peña,
% Javier Rodríguez-Sánchez, Alberto de Castro, Susana Marcos,
% Eduardo Martínez-Enríquez
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

switch origin
    case 'retina'
        for i = 1:Num_eyes
            [val(i) ind(i)] = max(dataZ_eigeneyes_pre(:,i));
            dataZ_eigeneyes_pre(:,i) = dataZ_eigeneyes_pre(:,i) - val(i);
        end

    case 'apex_cornea'
        for i = 1:Num_eyes
            [val(i) ind(i)] = min(dataZ_eigeneyes_pre(:,i));
            dataZ_eigeneyes_pre(:,i) = dataZ_eigeneyes_pre(:,i) - val(i);
        end

    case' '
end