addpath('../lib_matlab');

%
% Choosing nodes in parametric curve interpolation
% E T Y Lee
%
% computer-aided design
% volume 21 number 6 july/august 1989

PSET = {};

PSET{1} = [ 0, 0; 26, 24; 28, 24; 54, 0 ];
PSET{2} = [ 0, 0; 9, 39; 10, 40; 13, 40];
PSET{3} = [ 0, -0.15; 9.2, -0.15; 10, 0; 10, 0.5831];
PSET{4} = [ 0, -0.15; 9.2, -0.15; 10, 0; 10.5, 0.3];
PSET{5} = [ 0, -0.15; 9.2, -0.15; 10, 0; 9.5, 0.3];
PSET{6} = [ 0, 0; 10, 25; 10, 24; 11, 24.5; 33, 25];
PSET{7} = [ 0, 41; 10, 41; 31, 41; 40.8, 41; 41, 40.8; ...
            41, 31; 41, 10; 41, 0 ];
PSET{8}  = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0];
PSET{9}  = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31];
PSET{10} = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31; ...
             10, 41; 21, 41; 41, 21; 41, 10];
PSET{11} = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31; ...
             10, 41; 21, 41; 41, 21; 41, 10; 35, 41; 41, 35];
PSET{12} = [ 0, 0; 1.34, 5; 5, 8.66; 10, 10; 10.6, 10.4; 10.7, 12; ...
             10.7, 28.6; 10.8, 30.2; 11.4, 30.6; 19.6, 30.6; ...
             20.2, 30.2; 20.3, 28.6; 20.3, 12; 20.4, 10.4; ...
             21, 10; 26, 8.66; 29.66, 5; 31, 0 ];
PSET{13} = [ 0, 10; 2, 10; 3, 10; 5, 10; 6, 10; 8, 10; 9, 10.5; 11, 15; ...
             12, 50; 14, 60; 15, 85];
PSET{14} = [ 0, 1; 1, 1.1; 2, 1.1; 3, 1.2; 4, 1.3; 5, 7.2; 6, 3.1; ...
             7, 2.6; 8, 1.9; 9, 1.7; 10, 1.6];
PSET{15} = [ 2.5, 6.875; 5, 2.23; 10, 0.751; 15, 0.416; 20, 0.283; ...
             25, 0.219; 30, 0.182; 40, 0.143];
for k=1:15
  subplot( 3, 5, k );
%for k=1

  P = PSET{k};

  DP = (P(1:end-1,:)-P(2:end,:)).';
  % compute angles
  L     = sqrt(sum(DP.^2));
  theta = min(pi,acos( dot( DP(1:end-1), DP(2:end) )./(L(1:end-1).*L(2:end)) ));
  di    = L(2:end-1);
  dip1  = L(3:end);
  dim1  = L(1:end-2);
  thi   = theta(1:end-1);
  thip1 = theta(2:end);
  LL    = di.*( 1 + (3/2)*( thi.*dim1./(dim1+di) + thip1.*dip1./(di+dip1) ) );
  L1    = L(1).*( 1 + (3/2)*theta(1).*L(2)./(L(1)+L(2)) );
  LN    = L(end).*( 1 + (3/2)*theta(end).*L(end-1)./(L(end-1)+L(end)) );
  LL    = [L1,LL,LN];

  %T = sqrt(sum((P(1:end-1,:).'-P(2:end,:).').^2));
  %T = sqrt(T);
  T = [0,cumsum(LL)];
  T = T./T(end);

  if false
    S = SplineSet( 'cubic', T, P );
  else
    S = SplineVec();
    S.setup(P.');
    S.knots(linspace(0,1,length(T)));
    %S.knots(T);
    %S.chord();
    %S.centripetal();
    S.CatmullRom();
  end
  S.eval(0)
  S.eval(1)

  hold off;
  plot( P(:,1), P(:,2), ...
        '-o','Color','red', ...
        'MarkerSize',10,'MarkerFaceColor','blue','MarkerEdgeColor','green');
  hold on;
  tmin = S.tmin();
  tmax = S.tmax();
  t  = linspace(tmin,tmax,1000);
  PP = S.eval(t);
  plot( PP(1,:), PP(2,:), '-', 'Color', 'blue', 'Linewidth', 3);

  %axis equal
end
