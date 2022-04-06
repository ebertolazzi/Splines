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
%system('rmdir /S build');
if isunix
  system('cmake -Bbuild .');
  system('cd build; make -j 10; cd ..');
elseif ispc
  system('cmake -G "MinGW Makefiles" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -Bbuild .');
  system('cd build && make -j 10 && cd ..');
end
%system('cmake --build build --parallel 8');
cd(old_dir);

