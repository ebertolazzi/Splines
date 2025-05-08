%%-------------------------------------------------------------------------%%
%%                                                                         %%
%%  Copyright (C) 2020                                                     %%
%%                                                                         %%
%%         , __                 , __                                       %%
%%        /|/  \               /|/  \                                      %%
%%         | __/ _   ,_         | __/ _   ,_                               %%
%%         |   \|/  /  |  |   | |   \|/  /  |  |   |                       %%
%%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                      %%
%%                           /|                   /|                       %%
%%                           \|                   \|                       %%
%%                                                                         %%
%%      Enrico Bertolazzi                                                  %%
%%      Dipartimento di Ingegneria Industriale                             %%
%%      Università degli Studi di Trento                                   %%
%%      email: enrico.bertolazzi@unitn.it                                  %%
%%                                                                         %%
%%-------------------------------------------------------------------------%%

% Ottieni la dimensione dello schermo [left bottom width height]
scrn = get(0, 'ScreenSize');

% Calcola dimensioni all'80%
width  = round(scrn(3) * 0.8);
height = round(scrn(4) * 0.8);

% Centra la figura
left = round((scrn(3) - width) / 2);
bottom = round((scrn(4) - height) / 2);

% Crea la figura
figure('Position', [left bottom width height]);

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

if any(isnan(Z1(:)))
  fprintf('Found NaN on bicubic\n');
end
if any(isinf(Z1(:)))
  fprintf('Found Inf on bicubic\n');
end
if any(isnan(Z2(:)))
  fprintf('Found NaN on akima\n');
end
if any(isinf(Z2(:)))
  fprintf('Found Inf on akima\n');
end
if any(isnan(Z3(:)))
  fprintf('Found NaN on bilinear\n');
end
if any(isinf(Z3(:)))
  fprintf('Found Inf on bilinear\n');
end
if any(isnan(Z4(:)))
  fprintf('Found NaN on biquintic\n');
end
if any(isinf(Z4(:)))
  fprintf('Found Inf on biquintic\n');
end

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
