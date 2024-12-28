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

% Solicitar la iteraci�n deseada al usuario
iteracion_deseada = input('Introduce la iteraci�n que deseas representar: ');

% Filtrar los datos para la iteraci�n seleccionada
indices = (iteraciones == iteracion_deseada);
x = x_vals(indices);
E = E_vals(indices);
H = H_vals(indices);

% Graficar E(x)
figure(1);
clf;
plot(E, x, 'b-'); % Eje x: E, eje y: x
title(['Campo E para la iteraci�n ', num2str(iteracion_deseada)]);
xlabel('E');
ylabel('x');
grid on;

% Graficar H(x)
figure(2);
clf;
plot(H, x, 'r-'); % Eje x: H, eje y: x
title(['Campo H para la iteraci�n ', num2str(iteracion_deseada)]);
xlabel('H');
ylabel('x');
grid on;

%2,68, 134, 200

% Iteraciones que queremos representar
iteraciones_deseadas = [2, 68, 134, 200];

% Crear la figura con 4 subplots
figure;
clf;

for i = 1:length(iteraciones_deseadas)
    iteracion = iteraciones_deseadas(i);
    
    % Filtrar los datos para la iteraci�n actual
    indices = (iteraciones == iteracion);
    x = x_vals(indices);
    E = E_vals(indices);
    
    % Crear cada subplot
    subplot(2, 2, i);
    plot(E, x, 'b-', 'LineWidth', 1.5); % Eje x: E, eje y: x, con l�nea m�s gruesa
    title(['Campo E para la iteraci�n ', num2str(iteracion)], 'FontSize', 16);
    xlabel('E (V/m)', 'FontSize', 14);
    ylabel('y (m)', 'FontSize', 14);
    %legend(['Iteraci�n ', num2str(iteracion)], 'FontSize', 12, 'Location', 'best'); % A�adir leyenda
    grid on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 12); % Configuraci�n del grosor y tama�o de fuente del eje
end


