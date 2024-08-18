% Cp constants 
c1 = 0.22;
c2 = 116;
c3 = 0.4; 
c4 = 0;
c5 = 5; 
c6 = 12.5;

% Ranges for lambda and beta
lambda = linspace(0, 10, 200); 
beta_values = [0, 2, 5, 10, 20]; 


Cp = zeros(length(beta_values), length(lambda));

for i = 1:length(beta_values)
    beta = beta_values(i);
    lambda_i1 = 1 ./ (lambda + 0.08*beta) - 0.035 ./ (1 + beta^3);
    lambda_i = 1./lambda_i1;
    Cp(i,:) = c1 * (c2./lambda_i - c3*beta - c5) .* exp(-c6./lambda_i);
end

figure;
hold on;
colors = ['r', 'g', 'b', 'c', 'm']; 
for i = 1:length(beta_values)
    plot(lambda, Cp(i,:), colors(i), 'DisplayName', sprintf('\\beta = %d', beta_values(i)));
end
hold off;

xlabel('\lambda');
ylabel('C_p');
legend('show');
title('Performance Coefficient C_p as a Function of \lambda for Different \beta');
grid on;



%%
% Conform graficului, beta = 0:
lambda_opt = 6.33;
Cp_opt = 0.4382; % Coeficient de putere optime de 43.82%

% Conform graficului, pentru beta = 2:
lambda_opt2 = 7.33;
Cp_opt2 = 0.402; % Coeficient de putere optime de 40.02%

%%

ro = 1.225; %% air density (hg/m^3)
R = 5; %% blade radius (m)
A = pi*R^2; %% Turbine swept are (m^2)
wind_speed = 12; %% 12 m/s
Pm = 1/2 * Cp_opt * ro * A * wind_speed^3; 
% pt R = 5, Pm = 36.42 kW

omega_r = lambda_opt * wind_speed / R; 
% pt R = 5, wr = 15.19

%%


R = 5; 
rho = 1.225; 
A = pi * R^2; 

% Intervalul vitezei vantului
v = 0:0.5:25; % viteza vantului in m/s


P = 0.5 * rho * A * v.^3; % puterea in W


figure;
plot(v, P, 'LineWidth', 2);
xlabel('Viteza vantului (m/s)');
ylabel('Puterea (W)');
title('Relatia dintre viteza vantului si puterea produsa de o turbina eoliana cu R = 5 m');
grid on;

%%


theta_a = 0; % Phase A
theta_b = 120; % Phase B
theta_c = 240; % Phase C


theta_a_rad = deg2rad(theta_a);
theta_b_rad = deg2rad(theta_b);
theta_c_rad = deg2rad(theta_c);

V_magnitude = 1;


Vax = V_magnitude * cos(theta_a_rad);
Vay = V_magnitude * sin(theta_a_rad);

Vbx = V_magnitude * cos(theta_b_rad);
Vby = V_magnitude * sin(theta_b_rad);

Vcx = V_magnitude * cos(theta_c_rad);
Vcy = V_magnitude * sin(theta_c_rad);

% Plot the vectors
figure;
hold on;

% Plot thicker axes
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

% Add thicker lines for the x=0 and y=0 axes
plot([-1.5 1.5], [0 0], 'k', 'LineWidth', 1); % x=0 axis
plot([0 0], [-1.5 1.5], 'k', 'LineWidth', 1); % y=0 axis

quiver(0, 0, Vax, Vay, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5); % Phase A
quiver(0, 0, Vbx, Vby, 'g', 'LineWidth', 2.5, 'MaxHeadSize', 0.5); % Phase B
quiver(0, 0, Vcx, Vcy, 'b', 'LineWidth', 2.5, 'MaxHeadSize', 0.5); % Phase C

% Adding labels and grid
text(Vax, Vay, ' V_a', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(Vbx, Vby, ' V_b', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(Vcx, Vcy, ' V_c', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% Adding labels and grid
xlabel('X axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Vector Representation of V_{abc} with 120° Phase Shift', 'FontSize', 14);



% Plot arc between V_a and V_b
theta_arc = linspace(theta_a_rad, theta_b_rad, 100);
arc_x = 0.3 * cos(theta_arc);
arc_y = 0.3 * sin(theta_arc);
plot(arc_x, arc_y, 'k', 'LineWidth', 1);


% Set axes limits for better visibility
axis([-1.5 1.5 -1.5 1.5]);
grid on;
axis equal;
hold off;

%%

% Define phase angles
theta_a = 0; % Phase A
theta_b = 120; % Phase B
theta_c = 240; % Phase C

% Convert degrees to radians
theta_a_rad = deg2rad(theta_a);
theta_b_rad = deg2rad(theta_b);
theta_c_rad = deg2rad(theta_c);

% Define vector magnitudes (assuming equal magnitudes for simplicity)
V_magnitude = 1;

% Compute the x and y components of each vector
Vax = V_magnitude * cos(theta_a_rad);
Vay = V_magnitude * sin(theta_a_rad);

Vbx = V_magnitude * cos(theta_b_rad);
Vby = V_magnitude * sin(theta_b_rad);

Vcx = V_magnitude * cos(theta_c_rad);
Vcy = V_magnitude * sin(theta_c_rad);

% Clarke transformation
V_alpha = 2/3*(Vax-0.5*Vbx-0.5*Vcx);;
V_beta = 2/3*(sqrt(3)/2*Vby-sqrt(3)/2*Vcy);

% Plot the vectors in Clarke transformation
figure;
hold on;
% Plot thicker axes
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

% Add thicker lines for the \alpha=0 and \beta=0 axes
plot([-1.5 1.5], [0 0], 'k', 'LineWidth', 1); % \alpha=0 axis
plot([0 0], [-1.5 1.5], 'k', 'LineWidth', 1); % \beta=0 axis

quiver(0, 0, V_alpha, 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5); % Alpha axis
quiver(0, 0, 0, V_beta, 'g', 'LineWidth', 2.5, 'MaxHeadSize', 0.5); % Beta axis

% Adding labels and grid
text(V_alpha, 0, ' V_\alpha', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(0, V_beta, ' V_\beta', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% Adding labels and grid
xlabel('\alpha axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\beta axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Clarke Transformation Representation of V_{\alpha\beta}', 'FontSize', 14);

% Plot arc to show 90 degrees angle between alpha and beta axes
theta_arc = linspace(0, pi/2, 100);
arc_x = 0.3 * cos(theta_arc);
arc_y = 0.3 * sin(theta_arc);
plot(arc_x, arc_y, 'k', 'LineWidth', 1);
text(0.35 * cos(pi/4), 0.35 * sin(pi/4), '90°', 'FontSize', 12, 'FontWeight', 'bold');


axis([-1.5 1.5 -1.5 1.5]);
grid on;
axis equal;
hold off;


%%
% Define phase angles
theta_a = 0; % Phase A
theta_b = 120; % Phase B
theta_c = 240; % Phase C

% Convert degrees to radians
theta_a_rad = deg2rad(theta_a);
theta_b_rad = deg2rad(theta_b);
theta_c_rad = deg2rad(theta_c);

% Define vector magnitudes (assuming equal magnitudes for simplicity)
V_magnitude = 1;

% Compute the x and y components of each vector
Vax = V_magnitude * cos(theta_a_rad);
Vay = V_magnitude * sin(theta_a_rad);

Vbx = V_magnitude * cos(theta_b_rad);
Vby = V_magnitude * sin(theta_b_rad);

Vcx = V_magnitude * cos(theta_c_rad);
Vcy = V_magnitude * sin(theta_c_rad);

% Clarke transformation
V_alpha = Vax;
V_beta = (Vax + 2*Vbx) / sqrt(3);

% Define angular frequency and time
theta = pi / 6; % example angle of 30 degrees
cos_theta = cos(theta);
sin_theta = sin(theta);

% Rotate alpha and beta vectors
V_alpha_rotated = V_alpha * cos_theta - V_beta * sin_theta;
V_beta_rotated = V_alpha * sin_theta + V_beta * cos_theta;

% Plot the original vectors in Clarke transformation
figure;
hold on;

% Plot thicker axes
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

% Add thicker lines for the x=0 and y=0 axes
plot([-1.5 1.5], [0 0], 'k', 'LineWidth', 2); % x=0 axis
plot([0 0], [-1.5 1.5], 'k', 'LineWidth', 2); % y=0 axis

quiver(0, 0, V_alpha_rotated, V_beta_rotated, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Rotated alpha axis
quiver(0, 0, -V_beta_rotated, V_alpha_rotated, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Rotated beta axis

% Adding labels and grid
text(V_alpha_rotated, V_beta_rotated, ' V_{d}', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(-V_beta_rotated, V_alpha_rotated, ' V_{q}', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Adding labels and grid
xlabel('X axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Park Transformation Representation of V_{dq}', 'FontSize', 14);


% Plot arc to show theta angle between x-axis and V_alpha_rotated
theta_arc2 = linspace(0, theta, 100);
arc_x2 = 0.3 * cos(theta_arc2);
arc_y2 = 0.3 * sin(theta_arc2);
plot(arc_x2, arc_y2, 'k', 'LineWidth', 1, 'Color', 'm');

text(0.35 * cos(theta/2), 0.35 * sin(theta/2), '\theta', 'FontSize', 15, 'FontWeight', 'bold');


% Set axes limits for better visibility
axis([-1.5 1.5 -1.5 1.5]);
grid on;
axis equal;
hold off;


%%

% Cp constants 
c1 = 0.22;
c2 = 116;
c3 = 0.4; 
c4 = 0;
c5 = 5; 
c6 = 12.5;

% Ranges for lambda, beta, and wind speed
lambda = linspace(0, 10, 200); 
beta_values = [0, 2, 5, 10, 20]; 
wind_speeds = [6, 8, 10, 12, 14];  % Different wind speeds

Cp = zeros(length(beta_values), length(lambda));

% Calculate Cp for different beta values
for i = 1:length(beta_values)
    beta = beta_values(i);
    lambda_i1 = 1 ./ (lambda + 0.08*beta) - 0.035 ./ (1 + beta^3);
    lambda_i = 1./lambda_i1;
    Cp(i,:) = c1 * (c2./lambda_i - c3*beta - c5) .* exp(-c6./lambda_i);
end

% Parameters
ro = 1.225;  % Air density (kg/m^3)
R = 5;       % Blade radius (m)
A = pi * R^2; % Turbine swept area (m^2)

figure;
hold on;
colors = ['r', 'g', 'b', 'c', 'm']; 

% Calculate and plot power for different wind speeds
for j = 1:length(wind_speeds)
    wind_speed = wind_speeds(j);
    omega_r = lambda * wind_speed / R;
    Pm = 0.5 * ro * A * Cp(1,:) * wind_speed^3;
    plot(omega_r, Pm, colors(j), 'DisplayName', sprintf('v_w = %d m/s', wind_speeds(j)));
end

hold off;
xlabel('\omega_r (rad/s)');
ylabel('P_m (W)');
legend('show');
title('Mechanical Power P_m as a Function of Rotor Speed \omega_r for Different Wind Speeds');
grid on;

%%
% Cp constants 
c1 = 0.22;
c2 = 116;
c3 = 0.4; 
c4 = 0;
c5 = 5; 
c6 = 12.5;

% Ranges for wind speed and rotor speed
wind_speeds = linspace(4, 20, 200);  % Different wind speeds
omega_r_values = [5, 10, 15, 20, 25];  % Different rotor speeds

Cp = zeros(length(omega_r_values), length(wind_speeds));

% Parameters
ro = 1.225;  % Air density (kg/m^3)
R = 5;       % Blade radius (m)
A = pi * R^2; % Turbine swept area (m^2)
beta = 0;    % Fixed beta for simplicity

% Calculate Cp for different omega_r values
for i = 1:length(omega_r_values)
    omega_r = omega_r_values(i);
    lambda = omega_r * R ./ wind_speeds;
    lambda_i1 = 1 ./ (lambda + 0.08*beta) - 0.035 ./ (1 + beta^3);
    lambda_i = 1 ./ lambda_i1;
    Cp(i,:) = c1 * (c2 ./ lambda_i - c3 * beta - c5) .* exp(-c6 ./ lambda_i);
end

figure;
hold on;
colors = ['r', 'g', 'b', 'c', 'm']; 

% Plot Cp as a function of wind speed for different rotor speeds
for i = 1:length(omega_r_values)
    plot(wind_speeds, Cp(i,:), colors(i), 'DisplayName', sprintf('\\omega_r = %d rad/s', omega_r_values(i)));
end

hold off;
xlabel('Wind Speed (m/s)');
ylabel('C_p');
legend('show');
title('Performance Coefficient C_p as a Function of Wind Speed for Different Rotor Speeds');
grid on;

%%
% Parametri
wind_speed = 0:0.1:25; % Viteza vântului de la 0 la 25 m/s
cut_in_speed = 7;      % Cut-in speed (m/s)
power_limit_speed = 12; % Power limitation speed (m/s)
cut_out_speed = 20;    % Cut-out speed (m/s)

% Puterea nominală a turbinei (W)
P_nominal = 1500; % Exemplu: 1500 W

% Coeficient de putere, densitatea aerului, aria rotorului
Cp = 0.4; % Exemplu de coeficient de putere
ro = 1.225; % Densitatea aerului (kg/m^3)
R = 5; % Raza rotorului (m)
A = pi * R^2; % Aria rotorului

% Calcularea puterii produse la viteza de 12 m/s
P_limit = 0.5 * Cp * ro * A * power_limit_speed^3;

% Calcularea puterii produse
P_wind = zeros(size(wind_speed));

for i = 1:length(wind_speed)
    if wind_speed(i) < cut_in_speed
        P_wind(i) = 0;
    elseif wind_speed(i) >= cut_in_speed && wind_speed(i) < power_limit_speed
        % Folosim formula P = 0.5 * Cp * ro * A * v^3 pentru a calcula puterea
        P_wind(i) = 0.5 * Cp * ro * A * wind_speed(i)^3;
    elseif wind_speed(i) >= power_limit_speed && wind_speed(i) < cut_out_speed
        P_wind(i) = P_limit;
    else
        P_wind(i) = 0;
    end
end

% Generarea graficului
figure;
plot(wind_speed, P_wind, 'b', 'LineWidth', 3);
hold on;

% Adăugarea liniilor pentru vitezele de tăiere
y_limits = ylim;
plot([cut_in_speed cut_in_speed], y_limits, 'r--', 'LineWidth', 2);
plot([power_limit_speed power_limit_speed], y_limits, 'g--', 'LineWidth', 2);
plot([cut_out_speed cut_out_speed], y_limits, 'k--', 'LineWidth', 2);

% Adăugarea etichetelor
text(cut_in_speed, y_limits(2)*0.9, 'Cut-in Speed', 'Color', 'red', 'FontSize', 15, 'HorizontalAlignment', 'right');
text(power_limit_speed, y_limits(2)*0.9, 'Power Limitation', 'Color', 'green', 'FontSize', 15, 'HorizontalAlignment', 'right');
text(cut_out_speed, y_limits(2)*0.9, 'Cut-out Speed', 'Color', 'black', 'FontSize', 15, 'HorizontalAlignment', 'right');

% Ajustarea dimensiunii fontului pentru axe și titlu
xlabel('Viteza Vântului (m/s)', 'FontSize', 16);
ylabel('Puterea Produsă (W)', 'FontSize', 16);
title('Caracteristica Puterii Produse de Turbina Eoliană în Funcție de Viteza Vântului', 'FontSize', 17);
set(gca, 'FontSize', 12); % Ajustarea dimensiunii fontului pentru etichetele axelor
grid on;
hold off;

