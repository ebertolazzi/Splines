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
%%      Universita` degli Studi di Trento                                  %%
%%      email: enrico.bertolazzi@unitn.it                                  %%
%%                                                                         %%
%%-------------------------------------------------------------------------%%

X = 0:0.1:1;
Y = 0:0.25:2;

[XX,YY] = ndgrid(X,Y);

term1 = 0.75 * exp(-(9*XX-2).^2/4 - (9*YY-2).^2/4);
term2 = 0.75 * exp(-(9*XX+1).^2/49 - (9*YY+1)./10);
term3 = 0.5 * exp(-(9*XX-7).^2/4 - (9*YY-3).^2/4);
term4 = -0.2 * exp(-(9*XX-4).^2 - (9*YY-7).^2);

ZZ = term1 + term2 + term3 + term4;

ak = Spline2D('akima',X,Y,ZZ);

surf(XX,YY,ZZ);
view(145,-2);
set(gca,'Fontsize',16);

X = 0:0.01:1;
Y = 0:0.01:2;
[XX,YY] = ndgrid(X,Y);

ZZ = ak.eval(XX,YY);

surf(XX,YY,ZZ,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('akima');
