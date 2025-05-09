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

close all;

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

xx1 = [  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 ];
yy1 = [ 10, 11, 12, 11, 8, 12, 10.5, 15, 50, 60, 85 ];

fmt1 = {'LineWidth',3,'Color',[0.7,0.7,1]};
fmt2 = {'o','Color','red','MarkerSize',10,'MarkerFaceColor','yellow','MarkerEdgeColor','red'};
fmt3 = {'-','Color','black','Linewidth',2};
fmt4 = {'o','Color','blue','MarkerSize',10,'MarkerFaceColor','black','MarkerEdgeColor','black'};

subtype = 'pchip';

X = xx1;
Y = yy1;

%S = Spline1D('pchip',X,Y);
S = Spline1D('quintic',X,Y,subtype);

XX = X(1):(X(end)-X(1))/1000:X(end);
YY = S.eval(XX);

hold off;
plot( XX, YY, 'LineWidth', 2 );
hold on;

res = S.eval_y_min_max();
plot( res.x_min_pos, res.y_min, fmt2{:} );
plot( res.x_max_pos, res.y_max, fmt4{:} );

res = S.eval_y_min_max();
plot( res.x_min_pos, res.y_min, fmt2{:} );
plot( res.x_max_pos, res.y_max, fmt4{:} );
