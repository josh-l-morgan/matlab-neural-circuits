function []=print_pdf(file_name, width, height)

set(gcf, 'PaperPosition', [0.0 0.0 width height]);  % Centered left to right
set(gcf, 'PaperSize', [width height]);

set(gcf,'Renderer','painters'); % Seem to need these on linux to get vector pdf output
set(gcf,'RendererMode','manual');
print(gcf, '-dpdf', file_name );

