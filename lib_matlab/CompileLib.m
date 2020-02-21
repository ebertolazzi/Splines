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
  'SplinesUtils', ...
  'SplinesBivariate', ...
  'SplinesCinterface', ...
  'SplinesUnivariate'
};

LIB_SRCS = '';
LIB_OBJS = '';
for k=1:length(LIB_NAMES)
  LIB_SRCS = [ LIB_SRCS, ' ../src/', LIB_NAMES{k}, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES{k}, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES{k}, '.obj ' ];
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

  CMD = [ 'mex  -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -I../src -output ../lib_matlab/', N ];
  CMD = [ CMD, ' -largeArrayDims ../src_matlab_interface/mex_', N ];
  CMD = [ CMD, '.cc ', LIB_OBJS ];
  if isunix
    if ismac
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"'];
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"'];
    end
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end

for k=1:length(LIB_NAMES)
  if isunix
    delete([ LIB_NAMES{k}, '.o' ]);
  elseif ispc
    delete([ LIB_NAMES{k}, '.obj' ]);
  end
end

disp('----------------------- DONE ----------------------------');
