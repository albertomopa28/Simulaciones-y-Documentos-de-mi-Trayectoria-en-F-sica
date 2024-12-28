% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4�_Carrera\TFG\programas\VscodeFortran\Propagacion_unidireccional'
load salidas_unidireccional.txt;

filename = 'datos_unidireccional.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Leer el primer n�mero y asignarlo a la variable 'modo'
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo_modo = 'Pulso gaussiano';
elseif modo == 2
    titulo_modo = 'Excitaci�n +1 s�lo en l�nea 6: Modos TM';
else
    titulo_modo = 'Error';
end
fclose(fileID);

iter = salidas_unidireccional(:, 1);
t = salidas_unidireccional(:, 2);
E1 = salidas_unidireccional(:, 3);
H1 = salidas_unidireccional(:, 4);
E2 = salidas_unidireccional(:, 5);
H2 = salidas_unidireccional(:, 6);

n = length(t);
dt = t(2) - t(1);
c = 2.99792458e8;
dl = 2 * c * dt;
disp(dl);

hmin = 10 * dl;
fmaxn = c / hmin;
df = 1 / (n * dt);
f = (0:n-1) * df;
fmaxg = fmaxn;
nmf = round(fmaxg / df);

E1f = abs(fft(E1));
E2f = abs(fft(E2));
H1f = abs(fft(H1));
H2f = abs(fft(H2));


figure(1);
clf; % Nueva figura

% Primer subplot para la se�al E
subplot(2, 1, 1);
plot(iter, E1, 'b-'); 
hold on; % Permitir que otras curvas se a�adan
plot(iter, E2, 'g-');
hold off;
title('Campo E');
xlabel('Iteraci�n'); % Etiqueta del eje x
ylabel('E (V/m)'); % Etiqueta del eje y
legend('E1', 'E2'); % A�adir leyenda

% Segundo subplot para la se�al H
subplot(2, 1, 2);
plot(iter, H1, 'b-'); 
hold on; % Permitir que otras curvas se a�adan
plot(iter, H2, 'g-');
hold off;
title('Campo H');
xlabel('Iteraci�n'); % Etiqueta del eje x
ylabel('H (A/m)'); % Etiqueta del eje y
legend('H1', 'H2'); % A�adir leyenda

figure(2);
clf;

% Primer subplot para la se�al E
subplot(2, 1, 1);
hold on; % Permitir varias gr�ficas en la misma figura
plot(f(5:min(nmf, length(E1f)))*1e-6, E1f(5:min(nmf, length(E1f))));
plot(f(5:min(nmf, length(E2f)))*1e-6, E2f(5:min(nmf, length(E2f))));
hold off;
title('FFT Se�al E');
xlabel('Frecuencia (MHz)');
ylabel('mod(H)')
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ejes(3) = ampmin;
axis(ejes);

% Segundo subplot para la se�al H
subplot(2, 1, 2);
hold on; % Permitir varias gr�ficas en la misma figura
plot(f(5:min(nmf, length(H1f)))*1e-6, H1f(5:min(nmf, length(H1f))));
plot(f(5:min(nmf, length(H2f)))*1e-6, H2f(5:min(nmf, length(H2f))));
hold off;
title('FFT Se�al H');
xlabel('Frecuencia (MHz)');
ylabel('mod(H)')
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ejes(3) = ampmin;
axis(ejes);

