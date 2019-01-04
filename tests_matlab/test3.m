addpath('../lib_matlab');

X = -2:0.01:2;
Y = -2:0.01:2;
[XX,YY] = ndgrid(X,Y);
ZZ = peaks(XX,YY);
bc = Spline2D('bicubic',X,Y,ZZ);
ak = Spline2D('akima',X,Y,ZZ);
bl = Spline2D('bilinear',X,Y,ZZ);
bq = Spline2D('biquintic',X,Y,ZZ);

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);

X = -2:0.1:2;
Y = -2:0.1:2;
[XX,YY] = ndgrid(X,Y);
ZZ = peaks(XX,YY);

Z1 = bc.eval(XX,YY);
Z2 = ak.eval(XX,YY);
Z3 = bl.eval(XX,YY);
Z4 = bq.eval(XX,YY);

subplot(2,2,1);
surf(XX,YY,Z1,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bicubic');

subplot(2,2,2);
surf(XX,YY,Z2,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('akima');

subplot(2,2,3);
surf(XX,YY,Z3,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('bilinear');

subplot(2,2,4);
surf(XX,YY,Z4,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('biquintic');
