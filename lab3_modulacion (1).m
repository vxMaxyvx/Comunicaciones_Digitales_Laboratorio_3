%% =========================================================
%  LABORATORIO 3 – Modulación Pasabanda de Señales Binarias
%  CIT2111 – Universidad Diego Portales
%
%  Ejecutar sección por sección (Ctrl+Enter en cada %%)
%  o el script completo (F5).
%  Cada sección genera las figuras que necesitas para el informe.
% =========================================================
clear; clc; close all;

%% =========================================================
%  PARÁMETROS GLOBALES (modifica aquí para experimentar)
% =========================================================
R   = 100;          % Tasa de bits [bps]
Tb  = 1/R;          % Duración de bit [s]
Ac  = 1;            % Amplitud de portadora [V]
fc  = 1000;         % Frecuencia de portadora [Hz]  (debe cumplir fc >> R)
Df  = 200;          % Desviación FSK Δf [Hz]        (frecuencia de espacio = fc-Df, marca = fc+Df)
Nb  = 8;            % Número de bits a simular
fs  = 20*fc;        % Frecuencia de muestreo [Hz]  (20x portadora → buena resolución)
Ns  = fs*Tb;        % Muestras por bit

%% =========================================================
%  SECCIÓN 1 – SEÑAL MODULANTE m(t)
% =========================================================
% Generamos una secuencia binaria aleatoria y la repetimos
% Ns muestras por bit para construir la forma de onda continua.

rng(42);                            % Semilla fija → resultados reproducibles
bits = randi([0 1], 1, Nb);         % Secuencia binaria aleatoria (0 o 1)

% m(t) UNIPOLAR (para OOK): {0, 1}
m_unipolar = repelem(bits, Ns);

% m(t) POLAR (para FSK): {-1, +1}  →  2*bit - 1
m_polar = repelem(2*bits - 1, Ns);

t = (0 : length(m_unipolar)-1) / fs;   % Eje temporal

figure('Name','Señal modulante','NumberTitle','off','Position',[100 600 900 350]);
subplot(2,1,1)
  stairs(t*1e3, m_unipolar, 'b', 'LineWidth', 1.5);
  ylim([-0.2 1.4]); grid on;
  xlabel('Tiempo [ms]'); ylabel('Amplitud');
  title('Señal modulante m(t) UNIPOLAR – usada en OOK');
subplot(2,1,2)
  stairs(t*1e3, m_polar, 'r', 'LineWidth', 1.5);
  ylim([-1.4 1.4]); grid on;
  xlabel('Tiempo [ms]'); ylabel('Amplitud');
  title('Señal modulante m(t) POLAR – usada en FSK');

fprintf('=== PARÁMETROS ===\n');
fprintf('  Tasa de bits R  = %g bps\n', R);
fprintf('  Frecuencia fc   = %g Hz\n', fc);
fprintf('  Desviación Δf   = %g Hz (FSK)\n', Df);
fprintf('  Bits simulados  = %d\n', Nb);

%% =========================================================
%  SECCIÓN 2 – MODULACIÓN OOK (Actividad Previa)
% =========================================================
%
%  g(t) = Ac * m(t)          [envolvente compleja — real, unipolar]
%  s(t) = g(t)*cos(2π·fc·t)  = Ac * m(t) * cos(2π·fc·t)
%
%  Ancho de banda teórico (nulo a nulo):
%    BT_OOK = 2 * R = 2 / Tb

portadora = cos(2*pi*fc*t);         % cos(2π·fc·t)
g_ook     = Ac * m_unipolar;       % Envolvente compleja (real)
s_ook     = g_ook .* portadora;    % Señal modulada OOK

BT_OOK_teo = 2 * R;
fprintf('\n=== OOK ===\n');
fprintf('  Envolvente: g(t) = %.1f * m(t)\n', Ac);
fprintf('  BW teórico nulo-nulo: BT = 2R = %g Hz\n', BT_OOK_teo);

% --- Figura: dominio del tiempo ---
figure('Name','OOK – Dominio del tiempo','NumberTitle','off','Position',[100 200 900 500]);
subplot(3,1,1)
  stairs(t*1e3, m_unipolar, 'b', 'LineWidth',1.5); grid on;
  ylabel('m(t)'); title('Señal modulante unipolar m(t)');
subplot(3,1,2)
  plot(t*1e3, g_ook, 'k', 'LineWidth',1.2); grid on;
  ylabel('g(t)'); title('Envolvente compleja g(t) = A_c · m(t)');
subplot(3,1,3)
  plot(t*1e3, s_ook, 'r', 'LineWidth',0.8); grid on;
  xlabel('Tiempo [ms]'); ylabel('s(t)');
  title(['Señal OOK modulada  s(t) = g(t)·cos(2\pi·' num2str(fc) '·t)']);

%% =========================================================
%  SECCIÓN 3 – ESPECTRO DE LA SEÑAL OOK
% =========================================================
N    = length(s_ook);
f    = (-N/2 : N/2-1) * (fs/N);    % Eje de frecuencias centrado

% FFT de la envolvente compleja G(f)
G_ook = fftshift(fft(g_ook, N)) / fs;

% FFT de la señal modulada S(f)
S_ook = fftshift(fft(s_ook, N)) / fs;

figure('Name','OOK – Espectro','NumberTitle','off','Position',[100 100 950 500]);
subplot(2,1,1)
  plot(f, abs(G_ook), 'b', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|G(f)|');
  title('Espectro de la envolvente compleja |G(f)| – OOK');
  xlim([-3*R 3*R]);
  % Marcamos los nulos teóricos en ±R, ±2R, ...
  xline( R,'--r','LineWidth',1); xline(-R,'--r','LineWidth',1);
  legend('|G(f)|','Nulos en ±R');

subplot(2,1,2)
  plot(f, abs(S_ook), 'r', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|S(f)|');
  title(['Espectro señal OOK |S(f)| – BW teórico = 2R = ' num2str(BT_OOK_teo) ' Hz']);
  xlim([fc - 4*R, fc + 4*R]);
  xline(fc + R,'--k','LineWidth',1); xline(fc - R,'--k','LineWidth',1);
  legend('|S(f)|',['Nulos en f_c ± R = f_c ± ' num2str(R) ' Hz']);

% --- Medición del BW desde el espectro simulado ---
% Buscamos los nulos más cercanos a fc en la región positiva
f_pos    = f(f > 0);
S_pos    = abs(S_ook(f > 0));
[~, idx_fc] = min(abs(f_pos - fc));         % Índice de fc
% Buscar primer nulo a la derecha de fc
derecha  = S_pos(idx_fc:end);
[~, idx_nulo_der] = min(derecha(1:round(3*R*N/fs)));
f_nulo_der = f_pos(idx_fc + idx_nulo_der - 1);
BT_OOK_med = 2 * (f_nulo_der - fc);

fprintf('  BW medido (nulo-nulo desde espectro): %.1f Hz\n', BT_OOK_med);
fprintf('  Diferencia: %.2f %%\n', abs(BT_OOK_med - BT_OOK_teo)/BT_OOK_teo*100);

% --- TABLA DE RESULTADOS OOK (se imprime en consola) ---
fprintf('\n  ┌────────────────────────────────────┐\n');
fprintf('  │  TABLA COMPARATIVA – OOK           │\n');
fprintf('  ├──────────────┬──────────┬──────────┤\n');
fprintf('  │ Parámetro    │ Teórico  │ Medido   │\n');
fprintf('  ├──────────────┼──────────┼──────────┤\n');
fprintf('  │ fc (Hz)      │ %8g │   --     │\n', fc);
fprintf('  │ R  (bps)     │ %8g │   --     │\n', R);
fprintf('  │ BT (Hz)      │ %8g │ %8.1f │\n', BT_OOK_teo, BT_OOK_med);
fprintf('  └──────────────┴──────────┴──────────┘\n');

%% =========================================================
%  SECCIÓN 4 – MODULACIÓN FSK (Laboratorio Presencial)
% =========================================================
%
%  La FSK varía la frecuencia de la portadora según m(t):
%    f1 = fc - Df  → bit 0 (frecuencia de espacio)
%    f2 = fc + Df  → bit 1 (frecuencia de marca)
%
%  Señal modulada:
%    s(t) = Ac * cos(2π·fc·t + 2π·Δf · ∫m(σ)dσ)
%
%  Envolvente compleja:
%    g(t) = Ac · exp( j·2π·Δf · ∫m(σ)dσ )   → módulo constante Ac
%
%  Ancho de banda (regla de Carson):
%    BT_FSK = 2·(Δf + R)

% Integral de m(t) polar → fase acumulada
fase_acum = cumsum(m_polar) / fs;       % ∫m(σ)dσ  (integración numérica)

% Envolvente compleja (compleja)
g_fsk = Ac * exp(1j * 2*pi * Df * fase_acum);

% Señal modulada FSK (parte real de g(t)*e^{j·2π·fc·t})
s_fsk = real(g_fsk .* exp(1j * 2*pi * fc * t));

% Verificación alternativa directa (más intuitiva)
s_fsk2 = Ac * cos(2*pi*fc*t + 2*pi*Df*fase_acum);
% s_fsk y s_fsk2 deben ser idénticas (tolerancia numérica)

BT_FSK_teo = 2*(Df + R);
fprintf('\n=== FSK ===\n');
fprintf('  Envolvente: g(t) = %.1f · exp(j·2π·%g·∫m(σ)dσ)\n', Ac, Df);
fprintf('  |g(t)| = %.1f  (constante → modulación angular pura)\n', Ac);
fprintf('  BW teórico (Carson): BT = 2(Δf+R) = 2(%g+%g) = %g Hz\n', Df, R, BT_FSK_teo);

% --- Figura: dominio del tiempo ---
figure('Name','FSK – Dominio del tiempo','NumberTitle','off','Position',[150 200 900 600]);
subplot(4,1,1)
  stairs(t*1e3, m_polar, 'b', 'LineWidth',1.5); grid on; ylim([-1.4 1.4]);
  ylabel('m(t)'); title('Señal modulante polar m(t)  {-1,+1}');
subplot(4,1,2)
  plot(t*1e3, fase_acum, 'k', 'LineWidth',1.2); grid on;
  ylabel('\int m d\sigma'); title('Fase acumulada ∫m(σ)dσ');
subplot(4,1,3)
  plot(t*1e3, Ac*ones(size(t)), '--g', 'LineWidth',1);
  hold on;
  plot(t*1e3, real(g_fsk), 'm', 'LineWidth',1.2); grid on;
  ylabel('g(t)'); title('Envolvente compleja: Re{g(t)} (verde punteado = |g(t)| = constante)');
  legend('|g(t)| = A_c','Re{g(t)}');
subplot(4,1,4)
  plot(t*1e3, s_fsk, 'r', 'LineWidth',0.8); grid on;
  xlabel('Tiempo [ms]'); ylabel('s(t)');
  title(['Señal FSK modulada  (f_1=' num2str(fc-Df) ' Hz, f_2=' num2str(fc+Df) ' Hz)']);

%% =========================================================
%  SECCIÓN 5 – ESPECTRO DE LA SEÑAL FSK
% =========================================================
N2   = length(s_fsk);
f2   = (-N2/2 : N2/2-1) * (fs/N2);

% FFT envolvente compleja
G_fsk = fftshift(fft(g_fsk, N2)) / fs;

% FFT señal modulada
S_fsk = fftshift(fft(s_fsk, N2)) / fs;

figure('Name','FSK – Espectro','NumberTitle','off','Position',[150 100 950 500]);
subplot(2,1,1)
  plot(f2, abs(G_fsk), 'b', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|G(f)|');
  title('Espectro de la envolvente compleja |G(f)| – FSK');
  xlim([-3*(Df+R)  3*(Df+R)]);
  xline( Df,'--r','LineWidth',1,'Label','  +Δf');
  xline(-Df,'--r','LineWidth',1,'Label','-Δf  ');

subplot(2,1,2)
  plot(f2, abs(S_fsk), 'r', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|S(f)|');
  title(['Espectro señal FSK |S(f)| – BW Carson = ' num2str(BT_FSK_teo) ' Hz']);
  xlim([fc - 3*(Df+R), fc + 3*(Df+R)]);
  xline(fc + Df + R,'--k','LineWidth',1);
  xline(fc - Df - R,'--k','LineWidth',1);
  xline(fc + Df,    '--m','LineWidth',1,'Label','f_2 ');
  xline(fc - Df,    '--m','LineWidth',1,'Label','f_1 ');
  legend('|S(f)|','Límites BW Carson','','Frecuencias f_1, f_2');

% --- Medición BW desde espectro simulado ---
f_pos2   = f2(f2 > 0);
S_pos2   = abs(S_fsk(f2 > 0));
% Buscamos máximos locales cerca de fc-Df y fc+Df
[~, idx1] = min(abs(f_pos2 - (fc - Df)));
[~, idx2] = min(abs(f_pos2 - (fc + Df)));
f_left    = f_pos2(idx1);
f_right   = f_pos2(idx2);
BT_FSK_med = 2*(f_right - fc + R);   % Estimación por posición de lóbulos

% Buscar nulos externos al lóbulo derecho
lob_der = S_pos2(idx2 : min(idx2+round(2*R*N2/fs), end));
[~,idx_nulo] = min(lob_der(1:min(length(lob_der),round(R*N2/fs))));
f_nulo_ext = f_pos2(min(idx2 + idx_nulo - 1, length(f_pos2)));
BT_FSK_med2 = 2*(f_nulo_ext - fc);   % nulo a nulo

fprintf('  BW medido (posición de lóbulos): ~%.1f Hz\n', BT_FSK_med);
fprintf('  BW medido (nulo externo):        ~%.1f Hz\n', BT_FSK_med2);

fprintf('\n  ┌────────────────────────────────────────┐\n');
fprintf('  │  TABLA COMPARATIVA – FSK               │\n');
fprintf('  ├──────────────────┬──────────┬──────────┤\n');
fprintf('  │ Parámetro        │ Teórico  │ Medido   │\n');
fprintf('  ├──────────────────┼──────────┼──────────┤\n');
fprintf('  │ fc (Hz)          │ %8g │   --     │\n', fc);
fprintf('  │ Δf (Hz)          │ %8g │   --     │\n', Df);
fprintf('  │ R  (bps)         │ %8g │   --     │\n', R);
fprintf('  │ BT Carson (Hz)   │ %8g │ %8.1f │\n', BT_FSK_teo, BT_FSK_med);
fprintf('  └──────────────────┴──────────┴──────────┘\n');

%% =========================================================
%  SECCIÓN 6 – COMPARACIÓN OOK vs FSK
% =========================================================
figure('Name','Comparación OOK vs FSK','NumberTitle','off','Position',[200 150 1000 400]);
subplot(1,2,1)
  plot(f,  abs(S_ook), 'b', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|S(f)|');
  title('Espectro OOK'); xlim([fc-4*R, fc+4*R]);
  xline(fc+R,'--r'); xline(fc-R,'--r');

subplot(1,2,2)
  plot(f2, abs(S_fsk), 'r', 'LineWidth', 1.2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|S(f)|');
  title('Espectro FSK'); xlim([fc-3*(Df+R), fc+3*(Df+R)]);
  xline(fc+Df+R,'--k'); xline(fc-Df-R,'--k');
  xline(fc+Df,'--m'); xline(fc-Df,'--m');

fprintf('\n=== RESUMEN COMPARATIVO ===\n');
fprintf('  BW_OOK (teórico)  = 2·R          = 2·%g = %g Hz\n', R, BT_OOK_teo);
fprintf('  BW_FSK (teórico)  = 2·(Δf + R)   = 2·(%g+%g) = %g Hz\n', Df, R, BT_FSK_teo);
fprintf('  Relación BW_FSK/BW_OOK = %.2f\n', BT_FSK_teo/BT_OOK_teo);

%% =========================================================
%  SECCIÓN 7 – TRANSFORMADA DE FOURIER ANALÍTICA
%              (para graficar en el informe)
% =========================================================
% Esta sección grafica la forma analítica del espectro, que es
% la que debes mostrar en la parte teórica del informe.

f_base = linspace(-4*R, 4*R, 2000);   % eje de frecuencias normalizado

% |G_OOK(f)| ≈ Ac·Tb·|sinc(f·Tb)| + impulso en 0
% (ignoramos el impulso DC para la gráfica)
G_ook_analitico = abs(Ac * Tb * sinc(f_base * Tb));

% |G_FSK(f)| ≈ superposición de dos sinc desplazadas en ±Df
G_fsk_analitico = abs(0.5*Ac*Tb*sinc((f_base - Df)*Tb)) + ...
                  abs(0.5*Ac*Tb*sinc((f_base + Df)*Tb));

figure('Name','Espectros analíticos (para informe)','NumberTitle','off','Position',[250 100 1000 450]);
subplot(1,2,1)
  plot(f_base, G_ook_analitico, 'b', 'LineWidth', 2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|G(f)| [V/Hz]');
  title('|G(f)| OOK – Analítico');
  xline( R,'--r','Label','  R'); xline(-R,'--r','Label','-R');
  xline( 2*R,'--r'); xline(-2*R,'--r');

subplot(1,2,2)
  plot(f_base, G_fsk_analitico, 'r', 'LineWidth', 2); grid on;
  xlabel('Frecuencia [Hz]'); ylabel('|G(f)| [V/Hz]');
  title('|G(f)| FSK – Analítico');
  xline( Df,'--m','Label','  Δf'); xline(-Df,'--m','Label','-Δf');
  xline( Df+R,'--k'); xline(-Df-R,'--k');

fprintf('\n>> Todas las figuras generadas correctamente.\n');
fprintf('>> Guárdalas con: print(fig, ''nombre'', ''-dpng'', ''-r300'')\n');
fprintf('>> O usa: exportgraphics(gcf, ''nombre.png'', ''Resolution'', 300)\n');
