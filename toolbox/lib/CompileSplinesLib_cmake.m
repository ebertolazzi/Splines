clear all;
clear functions;
clear mex;
clc;

if ismac
  oldPath = getenv('PATH');
  newPath = strcat(oldPath, pathsep, '/usr/local/bin'); % on MAC
  setenv('PATH', newPath);
elseif isunix
elseif ispc
  oldPath = getenv('PATH');
  newPath = strcat(oldPath, pathsep, 'C:\Program Files\CMake\bin'); % on Windows
  setenv('PATH', newPath);
end

old_dir = cd(fileparts(which(mfilename)));
system('cmake -Bbuild .');
system('cmake --build build --parallel 8');
cd(old_dir);

