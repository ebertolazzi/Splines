clc;
clear functions;

NAMES = {
  'BaseHermiteWrapper', ...
  'SplineVecMexWrapper', ...
  'SplineSetMexWrapper', ...
  'Spline1DMexWrapper', ...
  'Spline2DMexWrapper' ...
};

LIB_NAMES = {
  'SplineAkima', ...
  'SplineAkima2D', ...
  'SplineBessel', ...
  'SplineBiCubic', ...
  'SplineBiQuintic', ...
  'SplineBilinear', ...
  'SplineConstant', ...
  'SplineCubic', ...
  'SplineCubicBase', ...
  'SplineHermite', ...
  'SplineLinear', ...
  'SplinePchip', ...
  'SplineQuintic', ...
  'SplineQuinticBase', ...
  'SplineSet', ...
  'SplineVec', ...
  'Splines', ...
  'Splines1D', ...
  'Splines2D', ...
  'SplinesUtils', ...
  'SplinesBivariate', ...
  'SplinesCinterface', ...
  'SplinesUnivariate'
};

LIB_SRCS = '';
LIB_OBJS = '';
for k=1:length(LIB_NAMES)
  n        = LIB_NAMES{k};
  LIB_SRCS = [ LIB_SRCS, ' ../src/', n, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, n, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, n, '.obj ' ];
  end
end
for k=1:length(NAMES)
  n        = NAMES{k};
  LIB_SRCS = [ LIB_SRCS, ' ../src_matlab_interface/mex_', n, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, n, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, n, '.obj ' ];
  end
end

[~,mexLoaded] = inmem('-completenames');


disp('---------------------------------------------------------');

CMD = 'mex -c  -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -largeArrayDims -I../src ';
if isunix
  if ismac
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
  else
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
  end
elseif ispc
end
CMD = [ CMD, LIB_SRCS ];

disp(CMD);
eval(CMD);

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex  -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -I../src -output ../lib_matlab/' ];
  CMD = [ CMD , N, ' -largeArrayDims ' ];
  if isunix
    CMD = [ CMD, N, '.o ', LIB_OBJS ];
  elseif ispc
    CMD = [ CMD, N, '.obj ', LIB_OBJS ];
  end
  if ismac
    CMD = [CMD, ' CXXFLAGS="\$CXXFLAGS -static-libgcc -static-libstdc++ -Wall -O2 -g"'];
  elseif isunix
    CMD = [CMD, ' CXXFLAGS="\$CXXFLAGS -static-libgcc -static-libstdc++ -Wall -O2 -g"'];
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end

for k=1:length(LIB_NAMES)
  n = LIB_NAMES{k};
  if isunix
    delete([ n, '.o' ]);
  elseif ispc
    delete([ n, '.obj' ]);
  end
end
for k=1:length(NAMES)
  n = NAMES{k};
  if isunix
    delete([ n, '.o' ]);
  elseif ispc
    delete([ n, '.obj' ]);
  end
end

disp('----------------------- DONE ----------------------------');
