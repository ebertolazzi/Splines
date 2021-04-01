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

xx0 = [  0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10 ]; 
yy0 = [ 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 ];

xx1 = [  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 ];
yy1 = [ 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 ];

xx2 = [  0,  2,  3,  5,  6,  8,    9, 11, 12, 14, 15 ];
yy2 = [ 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 ];

xx3 = [ 7.99, 8.09,       8.19,       8.7,      9.2,      10,       12,       15,       20       ];
yy3 = [ 0,    2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 ];

xx4 = [ 595,   635,   695,   795,   855,   875,   895,   915,   935,   985,   1035,  1075  ];
yy4 = [ 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 ];

xx5 = [ 0.11, 0.12, 0.15, 0.16 ];
yy5 = [ 0.0003, 0.0003, 0.0004, 0.0004 ];

XXX = { xx1, xx2, xx3, xx4, xx5, xx0 };
YYY = { yy1, yy2, yy3, yy4, yy5, yy0 };
LOC = {'northwest','northwest','southeast','northwest','northwest','northwest'};

close all;

subtype = 'pchip';

for k=1:6
  X = XXX{k};
  Y = YYY{k};

  pc = Spline1D('pchip',X,Y);
  qu = Spline1D('quintic',X,Y,subtype);

  XX = X(1):(X(end)-X(1))/1000:X(end);

  Y1 = pc.eval(XX);
  Y2 = qu.eval(XX);
	
  subplot(2,3,k);
  hold off;
  plot( X,  Y, 'o' );
  hold on;
  plot( XX, Y1, XX, Y2, 'LineWidth', 3 );

  legend('data','pchip','quintic');
  legend('boxoff');
  legend('Location',LOC{k});

end

figure();

for k=1:6
  X = XXX{k};
  Y = YYY{k};

  pc = Spline1D('pchip',X,Y);
  qu = Spline1D('quintic',X,Y,subtype);

  XX = X(1):(X(end)-X(1))/1000:X(end);

  Y1 = pc.eval_DD(XX);
  Y2 = qu.eval_DD(XX);
	
  subplot(2,3,k);
  plot( XX, Y1, XX, Y2, 'LineWidth', 3 );

  legend('pchip (y'''')','quintic (y'''')');
  legend('boxoff');
  legend('Location',LOC{k});

end
