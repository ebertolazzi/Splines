X = 0:0.1:1;
Y = 0:0.25:2;
[XX,YY] = ndgrid(X,Y);
ZZ = franke(XX,YY);

spline2d('ak','akima',X,Y,ZZ);

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);

X = 0:0.01:1;
Y = 0:0.01:2;
[XX,YY] = ndgrid(X,Y);

ZZ = spline2d('ak',X,Y);

surf(XX,YY,ZZ,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('akima');
