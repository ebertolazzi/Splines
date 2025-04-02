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
yy1 = [ 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 ];

xx2 = [  0,  2,  3,  5,  6,  8,    9, 11, 12, 14, 15 ];
yy2 = [ 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 ];

xx3 = [ 7.99, 8.09,       8.19,       8.7,      9.2,      10,       12,       15,       20       ];
yy3 = [ 0,    2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 ];

xx4 = [ 595,   635,   695,   795,   855,   875,   895,   915,   935,   985,   1035,  1075  ];
yy4 = [ 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 ];

XXX = { xx1, xx2, xx3, xx4 };
YYY = { yy1, yy2, yy3, yy4 };
LOC = {'northwest','northwest','southeast','northwest' };

fmt1 = {'LineWidth',3,'Color',[0.7,0.7,1]};
fmt2 = {'o','Color','red','MarkerSize',10,'MarkerFaceColor','yellow','MarkerEdgeColor','red'};
fmt3 = {'-','Color','black','Linewidth',2};
fmt4 = {'o','Color','blue','MarkerSize',10,'MarkerFaceColor','black','MarkerEdgeColor','black'};

close all;

subtype = 'pchip';

for k=1:4
  X = XXX{k};
  Y = YYY{k};

  pc = Spline1D('pchip',X,Y);
  qu = Spline1D('quintic',X,Y,subtype);

  XX = X(1):(X(end)-X(1))/1000:X(end);

  Y1 = pc.eval(XX);
  Y2 = qu.eval(XX);
	
  subplot(2,2,k);
  hold off;
  plot( X,  Y, 'o' );
  hold on;
  plot( XX, Y1, XX, Y2, 'LineWidth', 3 );

  legend('data','pchip','quintic');
  legend('boxoff');
  legend('Location',LOC{k});

  [ipos,xmi,ymi,jpos,xma,yma] = pc.eval_y_min_max();
  plot( [xmi,xma], [ymi,yma], fmt2{:} );

  [ipos,xmi,ymi,jpos,xma,yma] = qu.eval_y_min_max();
  plot( [xmi,xma], [ymi,yma], fmt4{:} );

end
