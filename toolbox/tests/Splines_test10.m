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

type = { 'bilinear', 'bicubic', 'biquintic', 'akima' };
set(gca,'Fontsize',16);

disableDefaultInteractivity(gca)

for k=1:4

  S = Spline2D(type{k});
  S.build('test_surf.json');

  x_min = S.x_min();
  x_max = S.x_max();
  y_min = S.y_min();
  y_max = S.y_max();

  X = x_min:0.05:x_max;
  Y = y_min:0.01:y_max;
  [XX,YY] = ndgrid(X,Y);

  ZZ = S.eval(XX,YY);

  ZZZ{k} = ZZ;

  if any(isnan(ZZ(:)))
    fprintf('Trovati NaN sulla superfice');
  end
  if any(isinf(ZZ(:)))
    fprintf('Trovati Inf sulla superfice');
  end

  subplot(2,2,k);

  surf(XX,YY,ZZ,'Linestyle',':');

  axis tight

  zlim([-1,30]);

  view(60,60);

  title(type{k});

end

