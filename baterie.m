%% Plot tensiunea bateriei în funcție de timp

plot(time, Vbat, "LineWidth", 1.5);
xlabel('Timp (s)');
ylabel('Tensiunea Bateriei (V)');
title('Tensiunea Bateriei în Funcție de Timp');
grid on;

% Plot curentul bateriei în funcție de timp
figure;
plot(time, current, "LineWidth", 1.5, "Color", 'm');
xlabel('Timp (s)');
ylabel('Curent (A)');
title('Curentul Bateriei în Funcție de Timp');
grid on;

%% Params Baterie

Cap_sys = 1111; % Ah
Cap_cell = 3.2; % Ah

V_nom_sys = 48;   % V
V_nom_cell = 3.6; % V

Nr_cells_series = round(V_nom_sys/V_nom_cell); % 13 celule serie
Nr_cells_parallel = round(Cap_sys/Cap_cell);   % 347 celule paralel

I_bat_impulse = 3.12; % A
T = 1020; % sec
Tau = 180; % sec

Rs_cell = (3.8-3.5)/I_bat_impulse;  % 0.0962 ohm
Rp_cell = (3.5-3.38)/I_bat_impulse; % 0.0385 ohm
Cp_cell = Tau/Rp_cell;              % 4772 Farad

Rs_sys = (Nr_cells_series * Rs_cell) / Nr_cells_parallel  % 0.0036 ohm
Rp_sys = (Nr_cells_series * Rp_cell) / Nr_cells_parallel  % 0.0014 ohm
Cp_sys = (Nr_cells_parallel * Cp_cell) / Nr_cells_series  % 124920 Farad

%% Params Second-order Model

Rp1 = 0.0099;
Cp1 = 1182;
Rp_sys1 = (Nr_cells_series * Rp1) / Nr_cells_parallel 
Cp_sys1 = (Nr_cells_parallel * Cp1) / Nr_cells_series 

Rp2 = 0.0375;
Cp2 = 3137;
Rp_sys2 = (Nr_cells_series * Rp2) / Nr_cells_parallel 
Cp_sys2 = (Nr_cells_parallel * Cp2) / Nr_cells_series  