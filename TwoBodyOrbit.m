%% TWO-BODY ORBIT PROPAGATOR (ELLIPTICAL OR CIRCULAR)
% Uses classical orbital elements -> Cartesian state
% Then numerically integrates r'' = -mu r / |r|^3

clear; clc; close all;

%% ---- Constants ----
mu = 3.986004418e14;      % Earth gravitational parameter [m^3/s^2]
Re = 6378.137e3;          % Earth radius [m]

%% ---- ORBIT DEFINITION (EDIT THESE) ----
% Perigee & apogee altitudes
hp = 100e4;               % perigee altitude [m]
ha = 5000e4;               % apogee altitude [m]

rp = Re + hp;
ra = Re + ha;

a = (rp + ra)/2;          % semi-major axis [m]
e = (ra - rp)/(ra + rp);  % eccentricity

% Orientation
i    = deg2rad(30);     % inclination [rad]
RAAN = deg2rad(30);       % right ascension of ascending node [rad]
argp = deg2rad(30);       % argument of perigee [rad]
nu   = deg2rad(0);        % true anomaly at start [rad]

%% ---- Initial Cartesian State ----
[r0, v0] = oe2rv(mu, a, e, i, RAAN, argp, nu);
x0 = [r0; v0];

fprintf('Initial radius   = %.3f km\n', norm(r0)/1e3);
fprintf('Initial speed    = %.3f km/s\n', norm(v0)/1e3);
fprintf('Eccentricity     = %.5f\n', e);

%% ---- Time Span ----
Torbit = 2*pi*sqrt(a^3/mu);    % orbital period [s]
tspan  = [0 2*Torbit];         % propagate 2 orbits

%% ---- Integrate Two-Body ODE ----
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t, X] = ode45(@(t,x) twobody_ode(t,x,mu), tspan, x0, opts);

r = X(:,1:3);    % position history [m]
v = X(:,4:6);    % velocity history [m/s]

%% =========================================================
%% OVERLAID STATE HISTORY PLOTS 
%% =========================================================

% Time in minutes
t_min = t / 60;

% ---- Acceleration history ----
a_hist = zeros(size(r));
for k = 1:length(t)
    a_hist(k,:) = -mu * r(k,:) / norm(r(k,:))^3;
end

% ---- Norms ----
r_norm = vecnorm(r,2,2);
v_norm = vecnorm(v,2,2);
a_norm = vecnorm(a_hist,2,2);

%% ===================== POSITION =====================
figure;
plot(t_min, r(:,1)/1e3, 'LineWidth',1.3); hold on
plot(t_min, r(:,2)/1e3, 'LineWidth',1.3)
plot(t_min, r(:,3)/1e3, 'LineWidth',1.3)
plot(t_min, r_norm/1e3, 'k', 'LineWidth',2.0)

grid on
xlabel('Time [min]')
ylabel('Position [km]')
title('Position Components and Magnitude vs Time')
legend('x','y','z','|r|','Location','best')

%% ===================== VELOCITY =====================
figure;
plot(t_min, v(:,1)/1e3, 'LineWidth',1.3); hold on
plot(t_min, v(:,2)/1e3, 'LineWidth',1.3)
plot(t_min, v(:,3)/1e3, 'LineWidth',1.3)
plot(t_min, v_norm/1e3, 'k', 'LineWidth',2.0)

grid on
xlabel('Time [min]')
ylabel('Velocity [km/s]')
title('Velocity Components and Magnitude vs Time')
legend('v_x','v_y','v_z','|v|','Location','best')

%% ===================== ACCELERATION =====================
figure;
plot(t_min, a_hist(:,1), 'LineWidth',1.3); hold on
plot(t_min, a_hist(:,2), 'LineWidth',1.3)
plot(t_min, a_hist(:,3), 'LineWidth',1.3)
plot(t_min, a_norm, 'k', 'LineWidth',2.0)

grid on
xlabel('Time [min]')
ylabel('Acceleration [m/s^2]')
title('Acceleration Components and Magnitude vs Time')
legend('a_x','a_y','a_z','|a|','Location','best')



%% ---- Energy Check (should be constant) ----
rmag = vecnorm(r,2,2);
vmag = vecnorm(v,2,2);
energy = 0.5*vmag.^2 - mu./rmag;

fprintf('Energy variation = %.3e J/kg\n', max(energy)-min(energy));

%% ---- Plot Orbit + Earth ----
figure;

% Draw Earth first (so orbit sits on top)
[xe, ye, ze] = sphere(80);  % resolution
hEarth = surf((Re*xe)/1e3, (Re*ye)/1e3, (Re*ze)/1e3);
set(hEarth, 'FaceAlpha', 0.25, 'EdgeColor', 'none');  % transparent, no mesh
hold on;

% Plot orbit
plot3(r(:,1)/1e3, r(:,2)/1e3, r(:,3)/1e3, 'LineWidth', 1.5)

% Axes / view
grid on; axis equal
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
title('Elliptical Two-Body Orbit (ECI)')

% Nice lighting (optional but looks good)
camlight headlight
lighting gouraud


%% =========================================================
%% FUNCTIONS
%% =========================================================

function xdot = twobody_ode(~, x, mu)
% Two-body point-mass dynamics
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    xdot = [v; a];
end

function [rI, vI] = oe2rv(mu, a, e, i, RAAN, argp, nu)
% Convert orbital elements to inertial Cartesian state

    p = a*(1 - e^2);   % semi-latus rectum

    % Perifocal position and velocity
    rPQW = (p/(1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    vPQW = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

    % Rotation matrix PQW -> ECI
    R = rot3(RAAN) * rot1(i) * rot3(argp);

    rI = R * rPQW;
    vI = R * vPQW;
end

function R = rot1(a)
% Rotation about x-axis
    R = [1 0 0;
         0 cos(a) -sin(a);
         0 sin(a)  cos(a)];
end

function R = rot3(a)
% Rotation about z-axis
    R = [ cos(a) -sin(a) 0;
          sin(a)  cos(a) 0;
          0        0     1];
end
