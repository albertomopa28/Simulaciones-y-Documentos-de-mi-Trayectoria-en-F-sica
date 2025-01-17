% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4�_Carrera\TFG\programas\VscodeFortran\Dispersion';
filename = 'video.txt';

% Leer el archivo de datos
fileID = fopen(filename, 'r');
data = textscan(fileID, '%f %f %f %f %f', 'HeaderLines', 1); % Ajustar si hay encabezado
fclose(fileID);

% Extraer los datos
iteraciones = data{1}; % N�mero de iteraci�n
tiempos = data{2};     % Tiempo asociado (en segundos)
x_vals = data{3};      % Posici�n en x
E_vals = data{4};      % Campo E
H_vals = data{5};      % Campo H

% Par�metros para el video
fps = 35;                          % Fotogramas por segundo
speed_factor = 8;                  % Factor de velocidad del video
frame_step = round(speed_factor);  % Paso entre iteraciones para el video
unique_iters = unique(iteraciones); % Iteraciones �nicas
unique_times = unique(tiempos); % Calcular los tiempos �nicos

% Crear el archivo de video
video_EH = VideoWriter('evolucion_EH_disp.avi');
video_EH.FrameRate = fps;
open(video_EH);

% Crear figura para graficar
fig = figure;

% Rango fijo para los ejes x (posici�n) y y (campos)
x_min = min(x_vals);
x_max = max(x_vals);
E_max = 1.1*max(E_vals); % Amplitud m�xima para E (basado en tu especificaci�n)
H_max = 1.1*max(H_vals); % Amplitud m�xima para H

% Generar el video iterando sobre los tiempos �nicos
for i = 1:frame_step:length(unique_iters)
    % Seleccionar datos para la iteraci�n actual
    iter_actual = unique_iters(i);
    indices = (iteraciones == iter_actual);
    
    x = x_vals(indices);
    E = E_vals(indices);
    H = H_vals(indices);

    % Verificar que los datos no est�n vac�os
    if isempty(x) || isempty(E) || isempty(H)
        continue; % Saltar si no hay datos
    end
    
    % Graficar E(x)
    subplot(2, 1, 1);
    plot(x, E, 'b');
    title(['Campo E en t = ', num2str(unique_times(i), '%.2e'), ' s, N=',num2str(iter_actual)]);
    xlabel('x m');
    ylabel('E V/m');
    xlim([x_min, x_max]);
    ylim([-E_max, E_max]);
    grid on;

    % Graficar H(x)
    subplot(2, 1, 2);
    plot(x, H, 'r');
    title(['Campo H en t = ', num2str(unique_times(i), '%.2e'), ' s, N=',num2str(iter_actual)]);
    xlabel('x m');
    ylabel('H A/m');
    xlim([x_min, x_max]);
    ylim([-H_max, H_max]);
    grid on;

    % Renderizar el gr�fico
    drawnow;

    % A�adir el fotograma al video
    frame = getframe(fig);
    writeVideo(video_EH, frame);
end

% Cerrar el archivo de video y la figura
close(video_EH);
close(fig);

disp('Video creado con �xito.');