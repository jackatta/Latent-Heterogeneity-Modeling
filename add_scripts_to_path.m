% Change the current folder to the folder of this m-file.
maindir = fileparts(which(mfilename));
cd(maindir);
pause(1);

% add functions to Matlab path if not already there
codedir = fullfile(maindir,'code');
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(codedir, pathCell));
else
  onPath = any(strcmp(codedir, pathCell));
end
if ~onPath
    addpath(genpath(codedir));
end

clear maindir codedir pathCell onPath