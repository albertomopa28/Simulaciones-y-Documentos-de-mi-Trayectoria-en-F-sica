% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4�_Carrera\TFG\programas\VscodeFortran\Guia\Guia_perfil_cuadrado';
filename = 'salidascuadrado.txt';

% Leer el archivo de datos
fileID = fopen(filename, 'r');
data = textscan(fileID, '%f %f %f %f %f', 'HeaderLines', 0);
fclose(fileID);

% Extraer los datos
iteraciones = data{1}; % N�mero de iteraci�n
tiempos = data{2};     % Tiempo asociado
x_vals = data{3};      % Valores de x
E_vals = data{4};      % Valores de E
H_vals = data{5};      % Valores de H

% Par�metros de tiempo y video
Delta_t = 4.8717039780534878E-12;  % Valor dado de ?t
fps = 25;                          % Fotogramas por segundo
speed_factor =1;                  % Factor de velocidad (mayor = m�s r�pido)
frame_step = round(speed_factor);  % Saltar iteraciones seg�n la velocidad
tiempos_unicos = unique(iteraciones) * Delta_t;  % Calcular los tiempos

% Crear el archivo de video
video_EH = VideoWriter('evolucion_EH.avi');
video_EH.FrameRate = fps;
open(video_EH);

% Obtener el rango global para y
y_vals_unicos = unique(x_vals);

% Crear la figura fuera del bucle y verificar si est� abierta
fig = figure(1);
for i = 1:frame_step:length(tiempos_unicos) % Ajustar paso del bucle
    % Verificar que la figura siga existiendo y est� abierta
    if ~isvalid(fig)
        disp('La figura se ha cerrado. Finalizando el proceso.');
        break;
    end
    
    % Seleccionar datos para la iteraci�n actual
    iteracion_actual = tiempos_unicos(i) / Delta_t;
    indices = (iteraciones == iteracion_actual);
    y = x_vals(indices);  % Usamos x como y
    E = E_vals(indices);
    H = H_vals(indices);
    
    % Verificar que no est�n vac�os los datos antes de graficar
    if isempty(y) || isempty(E) || isempty(H)
        continue;  % Saltar esta iteraci�n si no hay datos
    end
    
    % Graficar E(y)
    subplot(2,1,1);
    plot(E, y, 'b');
    title(['Campo E en t = ', num2str(tiempos_unicos(i), '%.2e'), ' s, N=', num2str(iteracion_actual)]);
    xlabel('E(y) V/m');
    ylabel('y m');
    
    % Ajustar los l�mites din�micamente en funci�n del valor m�ximo absoluto
    max_E = max(abs(E));
    if max_E > 0
        xlim([-1.1*max_E, 1.1*max_E]);
    else
        xlim([-1e-3, 1e-3]); % L�mite predeterminado en caso de max_E = 0
    end
    ylim([min(y_vals_unicos), max(y_vals_unicos)]);  % Rango fijo para y

    % Graficar H(y)
    subplot(2,1,2);
    plot(H, y, 'r');
    title(['Campo H en t = ', num2str(tiempos_unicos(i), '%.2e'), ' s, N=', num2str(iteracion_actual)]);
    xlabel('H(x) A/m');
    ylabel('y m');
    
    % Ajustar los l�mites din�micamente en funci�n del valor m�ximo absoluto
    max_H = max(abs(H));
    if max_H > 0
        xlim([-1.1*max_H, 1.1*max_H]);
    else
        xlim([-1e-5, 1e-5]); % L�mite predeterminado en caso de max_H = 0
    end
    ylim([min(y_vals_unicos), max(y_vals_unicos)]);  % Rango fijo para y

    % Asegurarse de que la figura se renderice
    drawnow;
    
    % A�adir el fotograma al video
    frame = getframe(fig);
    writeVideo(video_EH, frame);
end

% Cerrar el archivo de video y la figura
close(video_EH);
close(fig);

disp('Video creado con �xito.');



