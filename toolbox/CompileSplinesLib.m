clc;
clear functions;
[~,mexLoaded] = inmem('-completenames');

old_dir = cd(fileparts(which(mfilename)));

NAMES = {
  'SplineSetMexWrapper', ...
  'BaseHermiteWrapper', ...
  'SplineVecMexWrapper', ...
  'Spline1DMexWrapper', ...
  'Spline2DMexWrapper' ...
};

lst_cc = dir('./src/*.cc');

LIB_SRCS = '';
LIB_OBJS = '';

CMD = 'mex -c  -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -largeArrayDims -Isrc ';
if isunix
  if ismac
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -std=c++11 -Wall -O2 -g" '];
  else
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -std=c++11 -Wall -O2 -g" '];
  end
elseif ispc
end
CMD = [ CMD, LIB_SRCS ];

disp('---------------------------------------------------------');
for kk=1:length(lst_cc)
  name     = lst_cc(kk).name(1:end-3);
  LIB_SRCS = [ LIB_SRCS, ' src/', name, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, name, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, name, '.obj ' ];
  end
  CMD1 = [ CMD ' src/', name, '.cc' ];
  fprintf(1,'\n\nCompiling: %s.cc\n',name);
  disp(CMD1);
  eval(CMD1);
end

MROOT = matlabroot;

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'\n\nCompiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -Isrc -output bin/', N ];
  CMD = [ CMD, ' -largeArrayDims src_mex/mex_', N, '.cc ', LIB_OBJS ];
  if ismac
    CMD = [CMD, ' CXXFLAGS="\$CXXFLAGS -std=c++11 -Wall -O2 -g"'];
  elseif isunix
    % Workaround for MATLAB 2020 that force dynamic link with old libstdc++
    % solution: link with static libstdc++
    ARCH  = computer('arch');
    PATH1 = [MROOT, '/bin/', ARCH];
    PATH2 = [MROOT, '/extern/bin/', ARCH];
    CMD   = [ CMD, ...
      ' CXXFLAGS="\$CXXFLAGS -std=c++11 -Wall -O2 -g"' ...
      ' LDFLAGS="\$LDFLAGS -static-libgcc -static-libstdc++"' ...
      ' LINKLIBS="-L' PATH1 ' -L' PATH2 ' -lMatlabDataArray -lmx -lmex -lmat -lm "' ...
    ];
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end

for kk=1:length(lst_cc)
  name = lst_cc(kk).name(1:end-3);
  if isunix
    delete([ name, '.o' ]);
  elseif ispc
    delete([ name, '.obj' ]);
  end
end

cd(old_dir);

disp('----------------------- DONE ----------------------------');
