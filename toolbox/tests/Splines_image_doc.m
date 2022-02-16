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

P = [ 0, 0; 1.34, 5; 5, 8.66; 10, 10; 10.6, 10.4; 10.7, 12; ...
      10.7, 28.6; 10.8, 30.2; 11.4, 30.6; 19.6, 30.6; ...
      20.2, 30.2; 20.3, 28.6; 20.3, 12; 20.4, 10.4; ...
      21, 10; 26, 8.66; 29.66, 5; 31, 0 ];

S = SplineVec();
S.setup(P.');
S.centripetal();
S.CatmullRom();

% plot interpolation points
hold off;
plot( P(:,1), P(:,2), ...
      'o','Color','red', ...
      'MarkerSize',10,'MarkerFaceColor','blue','MarkerEdgeColor','green');
hold on;
tmin = S.tmin();
tmax = S.tmax();
t  = linspace(tmin,tmax,1000);
PP = S.eval(t);

% plot curve
plot( PP(1,:), PP(2,:), '-', 'Color', 'blue', 'Linewidth', 3);

save_png('exampleVec');

