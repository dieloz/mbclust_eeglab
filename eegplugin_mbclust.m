% eegplugin_mbclust() - EEGLAB plugin to cluster ICs based on finite 
% mixture model for probability densities
%
% Usage:
%   >> eegplugin_mbclust(fig, try_strings, catch_strings);
%
% Inputs:
%   fig           - [integer]  EEGLAB figure
%   try_strings   - [struct] "try" strings for menu callbacks.
%   catch_strings - [struct] "catch" strings for menu callbacks.
%
% See also:  pop_mbclust(), mbclust()
%
% Author: Diego Lozano-Soldevilla, Donders Institute, 04-Apr-2014 18:34:18

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Diego Lozano-Soldevilla, Donders Institute, dlozanosoldevilla@gmail.com
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers=eegplugin_mbclust( fig, try_strings, catch_strings);

 vers = 'mbc1.0';
    if nargin < 3
        error('eegplugin_mbclust requires 3 arguments');
    end
    
std = findobj(fig, 'tag', 'study');
uimenu(std, 'label', 'Cluster components by Model-Based Clustering', 'callback', ...
    [try_strings.check_ica '[STUDY ALLEEG]= pop_mbclust(STUDY,ALLEEG);' catch_strings.add_to_hist ], 'userdata', 'startup:off;study:on');
