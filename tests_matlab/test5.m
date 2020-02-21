addpath('../lib_matlab');

close all;

h = pi/16;
s = 0:h:4*pi;

c = 1/4;
d = 1/5;

X = [ (1+s*c).*cos(s) ];
Y = [ (1+s*c).*sin(s) ];
Z = [ s*d ];

fmt1 = {'LineWidth',3,'Color',[0.7,0.7,1]};
fmt2 = {'o','Color','red','MarkerSize',5,'MarkerFaceColor','yellow','MarkerEdgeColor','red'};
fmt3 = {'-','Color','black','Linewidth',2};
fmt4 = {'o','Color','blue','MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black'};

plot3(X,Y,Z,fmt2{:});
hold on

PNTS = [X(:).';Y(:).';Z(:).'];

npts = size(PNTS,2);

dp = sqrt(sum( (PNTS(:,1:end-1)-PNTS(:,2:end)).^2));
t  = [0,cumsum(dp)];
tt = t(1):(t(end)-t(1))/500000:t(end);


%type = 'cubic';
type = 'quintic';
%type = 'pchip';
%type = 'akima';
%type = 'bessel';
subtype = 'pchip';
%subtype = 'extrapolate';
Xs = Spline1D( type, t, PNTS(1,:), subtype );
Ys = Spline1D( type, t, PNTS(2,:), subtype );
Zs = Spline1D( type, t, PNTS(3,:), subtype );
plot3( Xs.eval(tt), Ys.eval(tt), Zs.eval(tt), fmt3{:} ) ;

figure();

plot(X,Y,fmt2{:});
hold on
plot( Xs.eval(tt), Ys.eval(tt), fmt3{:} );

figure();

dp = sqrt( (X(1:end-1)-X(2:end)).^2 +  (Y(1:end-1)-Y(2:end)).^2 + (Z(1:end-1)-Z(2:end)).^2 );
ss = [0,cumsum(dp)];

T_D   = [ Xs.eval_D(tt).'; Ys.eval_D(tt).'; Zs.eval_D(tt).' ];
T_DD  = [ Xs.eval_DD(tt).'; Ys.eval_DD(tt).'; Zs.eval_DD(tt).' ];
T_DDD = [ Xs.eval_DDD(tt).'; Ys.eval_DDD(tt).'; Zs.eval_DDD(tt).' ];

tmp = cross( T_D, T_DD );
kur = sqrt( sum( tmp.^2 ) ) ./  ( sum( T_D.^2 ) ).^1.5;
tor = dot( tmp, T_DDD ) ./ dot( tmp, tmp );

t1 = (s .^ 2);
t3 = (t1 + 2) .^ 2;
t4 = (c ^ 2);
t5 = (t4 .^ 2);
t13 = (d .^ 2);
t24 = sqrt((t5 .* t3 + t4 .* c .* (4 .* t1 .* s + 8 .* s) + t4 * (t1 .* (t13 + 6) + 4 .* t13 + 4) + c .* s .* (2 .* t13 + 4) + t13 + 1));
t28 = 2 * s .* c + t1 .* t4 + t13 + t4 + 1;
t29 = sqrt(t28);
kur1 = 0.1e1 ./ t29 ./ t28 .* t24;

t1 = (s .^ 2);
t3 = (c ^ 2);
t10 = (t1 + 2) .^ 2;
t11 = (t3 .^ 2);
t19 = (d ^ 2);
tor1 = 1 ./ (t11 .* t10 + t3 .* c .* (4 .* t1 .* s + 8 * s) + t3 .* (t1 .* (t19 + 6) + 4 * t19 + 4) + c .* s .* (2 * t19 + 4) + t19 + 1) .* (1 + t3 .* (t1 + 6) + 2 * s .* c) .* d;


subplot(2,1,1);
plot( tt, kur, ss, kur1 );
title('curvature');

subplot(2,1,2);
plot( tt, tor, ss, tor1 );
title('torsion');
ylim([-2,2]);

%plot( tt, Xs.eval_DD(tt) ) ;
%plot( tt, Xs.eval_D(tt) ) ;

figure();
%plot( tt, Xs.eval(tt), t, PNTS(1,:) );
plot( tt, Xs.eval_DD(tt) );
