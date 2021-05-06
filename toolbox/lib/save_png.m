function save_png( fname )
  % now our solution, export with two background colors
  exportgraphics(gcf, [ fname, '.pdf'], 'ContentType', 'image','ContentType','vector');
  system([ '/usr/local/bin/pdf2svg ', fname, '.pdf ', fname, '.svg'] );
  system([ '/usr/local/bin/pdftoppm -r 600 -png ', fname, '.pdf ', fname]);
  movefile([fname, '-1.png'], [ fname, '.png'])
end
