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
left   = round((scrn(3) - width) / 2);
bottom = round((scrn(4) - height) / 2);

% Crea la figura
figure('Position', [left bottom width height]);

X = [ 7.99, 8.09,       8.19,       8.7,      9.2,      10,       12,       15,       20       ];
Y = [ 0,    2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 ];

data.spline0.spline_type = 'cubic';
data.spline0.xdata       = X;
data.spline0.ydata       = Y;
data.spline1.spline_type = 'quintic';
data.spline1.spline_sub_type = 'pchip';
data.spline1.xdata       = X;
data.spline1.ydata       = Y;

spl = Spline1Dblend(data);

XX = X(1):(X(end)-X(1))/1000:X(end);
Y1 = spl.eval( XX, 0   );
Y2 = spl.eval( XX, 0.1 );
Y3 = spl.eval( XX, 0.3 );
Y4 = spl.eval( XX, 0.7 );
Y5 = spl.eval( XX, 0.9 );
Y6 = spl.eval( XX, 1   );

plot( X,  Y, 'o', ...
      XX, Y1, ...
      XX, Y2, ...
      XX, Y3, ...
      XX, Y4, ...
      XX, Y5, ...
      XX, Y6, 'LineWidth', 2 );
