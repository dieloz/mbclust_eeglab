% pop_mbclust() - select and run a clustering algorithm on components from an EEGLAB STUDY 
%               structure of EEG datasets. Clustering data should be prepared beforehand using 
%               pop_preclust() and/or std_preclust(). The number of clusters must be
%               specified in advance. If called in gui mode, the pop_clustedit() window
%               appears when the clustering is complete to display clustering results
%               and allow the user to review and edit them.
% Usage: 
%               >> STUDY = pop_mbclust( STUDY, ALLEEG); % pop up a graphic interface
%               >> STUDY = pop_mbclust( STUDY, ALLEEG, 'key1', 'val1', ...); % no pop-up
% Inputs:
%   STUDY       - an EEGLAB STUDY set containing some or all of the EEG sets in ALLEEG. 
%   ALLEEG      - a vector of loaded EEG dataset structures of all sets in the STUDY set.
%
% Optional Inputs:
%   'algorithm' - ['Model-Based Clustering'] cluster ICs based on finite mixture model for probability densities
%   'clus_num'  - [integer] the number of desired clusters (must be > 1) {default: 5}
%   'outliers'  - [integer] identify outliers further than the given number of standard
%                 deviations from any cluster centroid. CURRENTLY NOT AVAILABLE
%   'save'      - ['on' | 'off'] save the updated STUDY to disk {default: 'off'} 
%   'filename'  - [string] if save option is 'on', save the STUDY under this file name
%                    {default: current STUDY filename}
%   'filepath'  - [string] if save option is 'on', will save the STUDY in this directory 
%                    {default: current STUDY filepath}
% Outputs:
%   STUDY       - as input, but modified adding the clustering results.
%
% Graphic interface buttons:
%  "Clustering algorithm" - [list box] display/choose among the available clustering 
%                           algorithms. 
%  "Number of clusters to compute" - [edit box] the number of desired clusters (>1)
%
%  See also  pop_clustedit(), pop_preclust(), std_preclust(),
%
% Author: Diego Lozano-Soldevilla, Donders Institute, 04-Apr-2014 18:34:18
% 

% Copyright (C) Diego Lozano-Soldevilla, Donders Institute, dlozanosoldevilla@gmail.com
%
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

% Coding notes: Useful information on functions and global variables used.

function [STUDY, ALLEEG, command] = pop_mbclust(STUDY, ALLEEG, varargin)

command = '';
if nargin < 2
    help pop_mbclust;
    return;
end;

if isempty(STUDY.etc)
    error('No pre-clustering information, pre-cluster first!');
end;
if ~isfield(STUDY.etc, 'preclust')
    error('No pre-clustering information, pre-cluster first!');
end;
if isempty(STUDY.etc.preclust)
    error('No pre-clustering information, pre-cluster first!');
end;

if isempty(varargin) %GUI call
    
    % remove clusters below clustering level (done also after GUI)
    % --------------------------------------
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = [2:length(STUDY.cluster)];
    else
        for index = 2:length(STUDY.cluster)
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) & ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end;
        end;        
    end;
    
    if length(STUDY.cluster) > 2 & ~isempty(rmindex)
        resp = questdlg2('Clustering again will delete the last clustering results', 'Warning', 'Cancel', 'Ok', 'Ok');
        if strcmpi(resp, 'cancel'), return; end;
    end;
    
	alg_options = {'Model-Based Clustering'}; 
	set_outliers = ['set(findobj(''parent'', gcbf, ''tag'', ''outliers_std''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));'...
                            'set(findobj(''parent'', gcbf, ''tag'', ''std_txt''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));']; 
	algoptions = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''kmeans''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
                  
    strclust = '';
    if STUDY.etc.preclust.clustlevel > length(STUDY.cluster)
        STUDY.etc.preclust.clustlevel = 1;
    end;
    if STUDY.etc.preclust.clustlevel == 1
        strclust = [ 'Performing clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    else
        strclust = [ 'Performing sub-clustering on cluster ''' STUDY.cluster(STUDY.etc.preclust.clustlevel).name '''' ];
    end;
    
    numClustStr = '10';
    valalg = 1;
    
	clust_param = inputgui( { [1] [1] [1 1] [1 0.5 0.5 ]  }, ...
  { {'style' 'text'       'string' strclust 'fontweight' 'bold'  } {} ...
  {'style' 'text'       'string' 'Clustering algorithm:' 'fontweight' 'bold' } ...
  {'style' 'popupmenu'  'string' alg_options  'value' valalg 'tag' 'clust_algorithm'  'Callback' algoptions } ...
  {'style' 'text'       'string' 'Number of clusters to compute:' 'fontweight' 'bold' } ...
  {'style' 'edit'       'string' numClustStr 'tag' 'clust_num' } {}},...
  'pophelp(''pop_clust'')', 'Set clustering algorithm -- pop_clust()' , [] , 'normal', [ 1 .5 1 1 1]);
       
                          
	if ~isempty(clust_param)
        
        % removing previous cluster information
        % -------------------------------------
        if ~isempty(rmindex)
            fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
            STUDY.cluster(rmindex)          = [];
            STUDY.cluster(clustlevel).child = [];
            if clustlevel == 1 & length(STUDY.cluster) > 1
                STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
            end;
        end;
        
        clus_alg = alg_options{clust_param{1}};
        clus_num = str2num(clust_param{2});
                
        outliers = [];
        try
            clustdata = STUDY.etc.preclust.preclustdata;
        catch
            error('Error accesing preclustering data. Perform pre-clustering.');
        end;
        command = '[STUDY] = pop_mbclust(STUDY, ALLEEG,';
        
        if ~isempty(findstr(clus_alg, 'Model-Based Clustering')), clus_alg = 'Model-Based Clustering'; end;

        disp('Clustering ...');
        switch clus_alg
            case 'Model-Based Clustering'
                [bics,modelout,model,Z,IDX] = mbclust(clustdata, clus_num);
                [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm',  {'Model-Based Clustering', clus_num});
                command = sprintf('%s %s %d %s', command, '''algorithm'', ''Model-Based Clustering'',''clus_num'', ', clus_num, ',');
        end
        disp('Done.');
        
        % If save updated STUDY to disk
        save_on = 0; % old option to save STUDY
        if save_on
            command = sprintf('%s %s', command, '''save'', ''on'',');
            if ~isempty(clust_param{6})
                [filepath filename ext] = fileparts(clust_param{6});
                command = sprintf('%s%s%s%s%s%s', command, '''filename'', ''', [filename ext], ', ''filepath'', ''', filepath, ''');' );
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);
                STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [filename ext], 'filepath', filepath);
              else
                command(end:end+1) = ');';
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, command); 
                if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
                    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    STUDY = pop_savestudy(STUDY, ALLEEG);
                end
           end
       else
           command(end:end+1) = ');';
           STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);            
       end
           
       [STUDY com] = pop_clustedit(STUDY, ALLEEG); 
       command = [ command com];
	end
    
else %command line call
    % remove clusters below clustering level (done also after GUI)
    % --------------------------------------
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = [2:length(STUDY.cluster)];
    else
        for index = 2:length(STUDY.cluster)
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) & ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [ rmindex index ];
            end;
        end;        
    end;
    if ~isempty(rmindex)
        fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
        STUDY.cluster(rmindex)          = [];
        STUDY.cluster(clustlevel).child = [];
        if clustlevel == 1 & length(STUDY.cluster) > 1
            STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
        end;
    end;

    %default values
    clus_num = 20;
    save = 'off';
    filename = STUDY.filename;
    filepath = STUDY.filepath;
    outliers = Inf; % default std is Inf - no outliers
    
    if mod(length(varargin),2) ~= 0
        error('pop_clust(): input variables must be specified in pairs: keywords, values');
    end
    
    for k = 1:2:length(varargin)
        switch(varargin{k})
            case 'algorithm'
                algorithm = varargin{k+1};
            case 'clus_num'
                clus_num = varargin{k+1};
            case 'outliers'
                outliers =  varargin{k+1};                
            case 'save'
                save = varargin{k+1};
            case 'filename' 
                filename = varargin{k+1};
            case 'filepath'
                filepath = varargin{k+1};
        end
    end    
    if clus_num < 2
        clus_num = 2;
    end
    
    clustdata = STUDY.etc.preclust.preclustdata;
    switch lower(algorithm)
        case 'model-based clustering'
            [bics,modelout,model,Z,IDX] = mbclust(clustdata, clus_num);
            [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm',  {'Model-Based Clustering', clus_num});  
      otherwise
            disp('pop_clust: unknown algorithm return');
            return
    end
            % If save updated STUDY to disk
    if strcmpi(save,'on')
        if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
            STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
        else
            STUDY = pop_savestudy(STUDY);
        end
   end       
end
STUDY.saved = 'no';
