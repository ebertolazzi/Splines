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

V2 = zeros(10,10);
V2(2:5,2:5) = 3/7; % First step
V2(6:7,6:7) = [1/4 1/5; 1/5 1/4]; % Middle step
V2(8:9,8:9) = 1/2; % Last step
V2 = flip(V2,2);
X2 = 1:1:10;
Y2 = 1:1:10;

type = { 'bilinear', 'bicubic', 'biquintic', 'akima' };
set(gca,'Fontsize',16);

for k=1:4

  S = Spline2D(type{k});
  S.build(X2,Y2,V2);

  x_min = S.x_min();
  x_max = S.x_max();
  y_min = S.y_min();
  y_max = S.y_max();

  X = x_min:0.05:x_max;
  Y = y_min:0.01:y_max;
  [XX,YY] = ndgrid(X,Y);

  ZZ = S.eval(XX,YY);

  if any(isnan(ZZ(:)))
    fprintf('Trovati NaN sulla superfice');
  end
  if any(isinf(ZZ(:)))
    fprintf('Trovati Inf sulla superfice');
  end

  subplot(2,2,k);

  surf(XX,YY,ZZ,'Linestyle',':');
  title(type{k});

  axis tight

  zlim([-0.5,1]);

  view(60,60);

end

