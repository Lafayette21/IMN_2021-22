matrix_g1 = readmatrix('V16.txt');
mapa = pcolor(matrix_g1');
set(mapa, 'EdgeColor', 'none');
colormap turbo
colorbar
title('k=16','FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlabel('x','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold') 
ylabel('y','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold')

figure;

matrix_g2 = readmatrix('V8.txt');
mapa = pcolor(matrix_g2');
set(mapa, 'EdgeColor', 'none');
colormap turbo
colorbar
title('k=8','FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlabel('x','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold') 
ylabel('y','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold')

figure;

matrix_g3 = readmatrix('V4.txt');
mapa = pcolor(matrix_g3');
set(mapa, 'EdgeColor', 'none');
colormap turbo
colorbar
title('k=4','FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlabel('x','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold') 
ylabel('y','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold')

figure;

matrix_g4 = readmatrix('V2.txt');
mapa = pcolor(matrix_g4');
set(mapa, 'EdgeColor', 'none');
colormap turbo
colorbar
title('k=2','FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlabel('x','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold') 
ylabel('y','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold')

figure;

matrix_g5 = readmatrix('V1.txt');
mapa = pcolor(matrix_g5');
set(mapa, 'EdgeColor', 'none');
colormap turbo
colorbar
title('k=1','FontSize', 20, 'Color', 'k', 'FontWeight', 'bold');
xlabel('x','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold') 
ylabel('y','FontSize', 15, 'Color', 'k', 'FontWeight', 'bold')