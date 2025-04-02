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

xx1 = [  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 ];
%yy1 = [  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 ];
yy1 = [ 10, 11, 12, 11, 8, 12, 10.5, 15, 50, 60, 85 ];

close all;

subtype = 'pchip';

X = xx1;
Y = yy1;

S = Spline1D('pchip',X,Y);
%S = Spline1D('quintic',X,Y,subtype);

[coeffs,nodes] = S.eval_coeffs();

coeffs

nodes