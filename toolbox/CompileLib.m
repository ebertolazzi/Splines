clc;
clear functions;
[~,mexLoaded] = inmem('-completenames');

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
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
  else
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
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
  fprintf(1,'Compiling: %s.cc\n',name);
  eval(CMD1);
end

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex  -DSPLINES_DO_NOT_USE_GENERIC_CONTAINER -Isrc -output bin/', N ];
  CMD = [ CMD, ' -largeArrayDims src_mex/mex_', N ];
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

for kk=1:length(lst_cc)
  name = lst_cc(kk).name(1:end-3);
  if isunix
    delete([ name, '.o' ]);
  elseif ispc
    delete([ name, '.obj' ]);
  end
end

disp('----------------------- DONE ----------------------------');
