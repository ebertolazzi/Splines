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
ZZ = franke(XX,YY);

ak = Spline2D('akima',X,Y,ZZ);

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);

X = 0:0.01:1;
Y = 0:0.01:2;
[XX,YY] = ndgrid(X,Y);

ZZ = ak.eval(XX,YY);

surf(XX,YY,ZZ,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);
title('akima');
