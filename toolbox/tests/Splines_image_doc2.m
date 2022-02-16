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

X = -2:0.01:2;
Y = -2:0.01:2;
[XX,YY] = ndgrid(X,Y);
ZZ  = peaks(XX,YY);
spl = Spline2D('bicubic',X,Y,ZZ);

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);

X = -2:0.1:2;
Y = -2:0.1:2;
[XX,YY] = ndgrid(X,Y);

ZZ = spl.eval(XX,YY);

surf(XX,YY,ZZ,'Linestyle',':'), view(145,40), set(gca,'Fontsize',16);

save_png('exampleSurf');