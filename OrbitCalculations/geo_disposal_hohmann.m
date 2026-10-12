%% GEO End-of-Mission Disposal: Hohmann Transfer to Graveyard Orbit
% Solves the two-burn Hohmann transfer (first burn at periapsis of the
% initial orbit) and draws an interactive 3D presentation figure.
clear; clc; close all;

%% 1) Constants and orbit data
mu   = 398600.4418;   % Earth gravitational parameter [km^3/s^2]
Re   = 6378.137;      % Earth equatorial radius [km]
rGEO = 42164.17;      % geostationary radius, shown as a reference ring [km]

a1  = 42464;  e1 = 0.005;   % initial orbit
a2  = 42564;                % final orbit (circular, e = 0)
inc = 1;                    % inclination of all orbits [deg]

%% 2) Hohmann transfer solution
rp1 = a1*(1 - e1);                  % initial periapsis radius (burn 1) [km]
ra1 = a1*(1 + e1);                  % initial apoapsis radius [km]
aT  = (rp1 + a2)/2;                 % transfer semi-major axis [km]
eT  = (a2 - rp1)/(a2 + rp1);        % transfer eccentricity

vis = @(r,a) sqrt(mu*(2./r - 1./a));    % vis-viva speed [km/s]
dV1 = vis(rp1,aT) - vis(rp1,a1);        % burn 1 at initial periapsis [km/s]
dV2 = vis(a2,a2)  - vis(a2,aT);         % burn 2 at transfer apoapsis [km/s]
dVtot = abs(dV1) + abs(dV2);            % total delta-v [km/s]
TOF   = pi*sqrt(aT^3/mu);               % half-period transfer time [s]

fprintf('Burn 1 (periapsis):  dV1 = %+.3f m/s\n', dV1*1e3);
fprintf('Burn 2 (apoapsis):   dV2 = %+.3f m/s\n', dV2*1e3);
fprintf('Total delta-v:       %.3f m/s\n', dVtot*1e3);
fprintf('Transfer time:       %.2f hours\n', TOF/3600);

%% 3) Orbit geometry for plotting
% The orbits are only a few hundred km apart at ~42,000 km, so each apsis
% radius is pushed away from the GEO ring by a factor k. Every orbit is then
% drawn as an exact Keplerian ellipse (Earth at the focus) through its
% exaggerated periapsis and apoapsis, so the shapes, apsis locations, and
% tangency at both burns are all kept. With k = 1 these are the true orbits.
% RAAN and argument of periapsis were not given, so both are set to 0
% (periapsis on the +X axis). They only rotate the picture, not the dV.
k = 25;                                          % set k = 1 for true scale
scaleR = @(r) rGEO + k*(r - rGEO);               % exaggerate distance from GEO

rpP = scaleR(rp1);  raP = scaleR(ra1);  rfP = scaleR(a2);   % plotted radii
aP1 = (rpP + raP)/2;  eP1 = (raP - rpP)/(raP + rpP);        % plotted initial orbit
aPT = (rpP + rfP)/2;  ePT = (rfP - rpP)/(rfP + rpP);        % plotted transfer orbit

Rx      = [1 0 0; 0 cosd(inc) -sind(inc); 0 sind(inc) cosd(inc)];  % tilt by inclination
radius  = @(a,e,nu) a*(1 - e^2)./(1 + e*cos(nu));                   % conic equation
toXYZ   = @(r,nu) Rx*[r.*cos(nu); r.*sin(nu); zeros(size(nu))];    % orbit plane -> inertial
onOrbit = @(a,e,deg) toXYZ(radius(a,e,deg*pi/180), deg*pi/180);    % point at true anomaly [deg]

nu  = linspace(0, 2*pi, 721);                    % full orbit, true anomaly [rad]
nuT = linspace(0, pi, 361);                      % transfer = half ellipse (periapsis -> apoapsis)
P1  = toXYZ(radius(aP1,eP1,nu), nu);             % initial orbit
PT  = toXYZ(radius(aPT,ePT,nuT), nuT);           % transfer orbit
P2  = toXYZ(rfP*ones(size(nu)), nu);             % final (graveyard) orbit
PG  = toXYZ(rGEO*ones(size(nu)), nu);            % GEO reference ring
B   = [toXYZ(rpP, 0), toXYZ(rfP, pi)];           % burn 1 and burn 2 locations

%% 4) Figure with a starfield background
cInit  = [0.35 0.62 1.00];   % initial orbit  (blue)
cXfer  = [1.00 0.55 0.18];   % transfer orbit (orange)
cFinal = [0.25 0.88 0.62];   % final orbit    (green)
cBurn  = [1.00 0.92 0.45];   % burn delta-v arrows (yellow)

fig = figure('Color','k', 'Position',[100 100 1280 800]);

starAx = axes('Position',[0 0 1 1]);             % flat background layer for stars
rng(7);  nStars = 900;                           % fixed seed = same stars every run
scatter(starAx, rand(nStars,1), rand(nStars,1), 0.5 + 6*rand(nStars,1).^6, ...
        (0.35 + 0.65*rand(nStars,1))*[1 1 1], 'filled');
axis(starAx, [0 1 0 1]);  axis(starAx, 'off');

ax = axes('Position',[0 0.02 1 0.86], 'Color','none');   % main 3D scene on top
hold(ax, 'on');  axis(ax, 'equal', 'off');

%% 5) Earth with NASA Blue Marble texture and a thin atmosphere
texFile = 'earth_texture.jpg';
texURL  = 'https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_2048.jpg';
try
    if ~isfile(texFile), websave(texFile, texURL); end    % downloads once, then reuses the file
    earthImg = flipud(imread(texFile));                   % flip so north is up on the sphere
catch
    earthImg = repmat(reshape([0.10 0.30 0.60],1,1,3), 2, 2);   % offline: plain blue Earth
end

[xs,ys,zs] = sphere(150);
surf(Re*xs, Re*ys, Re*zs, 'CData',earthImg, 'FaceColor','texturemap', 'EdgeColor','none', ...
     'AmbientStrength',0.25, 'DiffuseStrength',1, 'SpecularStrength',0.3, 'SpecularExponent',15);
surf(1.08*Re*xs, 1.08*Re*ys, 1.08*Re*zs, 'FaceColor',[0.35 0.60 1.00], 'EdgeColor','none', ...
     'FaceAlpha',0.15, 'AmbientStrength',0.2);                          % atmosphere haze

%% 6) Orbits (with a soft glow)
% Glow = a wide, faint copy of a line drawn underneath it (NaN keeps the patch open)
glow = @(P,c,w,a) patch([P(1,:) NaN], [P(2,:) NaN], [P(3,:) NaN], 'k', ...
                        'FaceColor','none', 'EdgeColor',c, 'EdgeAlpha',a, 'LineWidth',w);

hG = plot3(PG(1,:), PG(2,:), PG(3,:), ':', 'Color',[0.55 0.55 0.60], 'LineWidth',1);
h1 = plot3(P1(1,:), P1(2,:), P1(3,:), '--', 'Color',cInit, 'LineWidth',1.5);

glow(P2, cFinal, 10, 0.18);
h2 = plot3(P2(1,:), P2(2,:), P2(3,:), '-', 'Color',cFinal, 'LineWidth',1.8);

glow(PT, cXfer, 18, 0.14);  glow(PT, cXfer, 9, 0.30);
hT = plot3(PT(1,:), PT(2,:), PT(3,:), '-', 'Color',cXfer, 'LineWidth',2.8);

%% 7) Direction-of-motion arrows (all orbits are prograde: counterclockwise seen from north)
% Each arrowhead sits on its orbit and points along the orbit's tangent
headLen = 3600;                                  % arrowhead length [plot km]
arrowOn = @(a,e,deg,c) drawCone(onOrbit(a,e,deg) + headLen/2*unitVec(onOrbit(a,e,deg+1) - onOrbit(a,e,deg-1)), ...
                                onOrbit(a,e,deg+1) - onOrbit(a,e,deg-1), headLen, c);
arrowOn(aP1, eP1, 270, cInit);
arrowOn(aPT, ePT,  45, cXfer);   arrowOn(aPT, ePT, 150, cXfer);
arrowOn(rfP,   0,  60, cFinal);  arrowOn(rfP,   0, 300, cFinal);

%% 8) Burn points and delta-v arrows
% At an apsis the velocity is perpendicular to the radius, so the thrust
% direction is +/- the flight direction. Arrow length is proportional to |dV|.
vHat    = @(nu) Rx*[-sin(nu); cos(nu); 0];               % flight direction at an apsis
burnDir = [sign(dV1)*vHat(0), sign(dV2)*vHat(pi)];      % burn 1 retrograde, burn 2 prograde
burnLen = 4000*abs([dV1 dV2])*1e3;                       % 4000 plot km per m/s of dV

for b = 1:2
    tip = B(:,b) + burnLen(b)*burnDir(:,b);
    hDV = plot3([B(1,b) tip(1)], [B(2,b) tip(2)], [B(3,b) tip(3)], 'Color',cBurn, 'LineWidth',3);
    drawCone(tip, burnDir(:,b), headLen, cBurn);
end

scatter3(B(1,:), B(2,:), B(3,:), 700, cXfer, 'filled', 'MarkerFaceAlpha',0.25);   % burn glow
scatter3(B(1,:), B(2,:), B(3,:), 60, 'w', 'filled');                               % burn core

% Burn labels: anchored just outside the outermost orbit on each side of the
% line of apsides, so they never sit on top of an orbit line
gray = '\color[rgb]{0.65 0.70 0.80}';
L1 = toXYZ(rfP + 2500, 0);                       % beyond the graveyard orbit (burn 1 side)
L2 = toXYZ(raP + 2500, pi);                      % beyond the initial apoapsis (burn 2 side)
text(L1(1), L1(2), L1(3), {sprintf('\\bf\\DeltaV_1 = %+.2f m/s', dV1*1e3), [gray '\rmretrograde burn']}, ...
     'Color','w', 'FontSize',14, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
text(L2(1), L2(2), L2(3), {sprintf('\\bf\\DeltaV_2 = %+.2f m/s', dV2*1e3), [gray '\rmprograde burn']}, ...
     'Color','w', 'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','middle');

%% 9) Simple satellite model on the transfer arc
nuS = 115*pi/180;                                 % where to draw it [rad]
pS  = toXYZ(radius(aPT,ePT,nuS), nuS);            % satellite position
uR  = pS/norm(pS);                                % radial direction
uN  = Rx*[0;0;1];                                 % orbit normal
uT  = cross(uN, uR);                              % along-track direction
M   = [uR uT uN];                                 % local -> inertial rotation

s = 1800;                                         % icon size (not to scale) [km]
busV   = s*[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
busF   = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
panelV = s*[0 -0.8 1.3; 0 0.8 1.3; 0 0.8 4.5; 0 -0.8 4.5];   % one solar array (north-south, like a GEO sat)

toWorld = @(V) (M*V' + pS)';                      % move local vertices into place
patch('Vertices',toWorld(busV), 'Faces',busF, 'FaceColor',[0.90 0.70 0.30], ...
      'EdgeColor',[0.3 0.2 0.1], 'AmbientStrength',0.6, 'SpecularStrength',0.9);
for side = [1 -1]                                 % two arrays, north and south
    patch('Vertices',toWorld(panelV.*[1 1 side]), 'Faces',1:4, 'FaceColor',[0.10 0.20 0.50], ...
          'EdgeColor',[0.55 0.65 0.90], 'AmbientStrength',0.6, 'SpecularStrength',0.8);
end

%% 10) Lighting, camera, and mouse rotation
light(ax, 'Position',[0.3 -1 0.7]);              % "Sun" direction (fixed in space)
lighting(ax, 'gouraud');
axis(ax, 'vis3d');                               % keep the scale fixed while rotating
view(ax, -15, 50);  camzoom(ax, 1.5);            % starting view angle and zoom

r = rotate3d(fig);  r.Enable = 'on';             % click and drag to rotate the scene
setAllowAxesRotate(r, starAx, false);            % keep the star background still

%% 11) Legend, title, and footnote (fixed on screen)
legend([h1 hT h2 hG hDV], { ...
    sprintf('Initial orbit  (a = %.0f km, e = %.3f)', a1, e1), ...
    sprintf('Hohmann transfer  (a = %.0f km, e = %.4f, %.2f h)', aT, eT, TOF/3600), ...
    sprintf('Graveyard orbit  (a = %.0f km, circular)', a2), ...
    'GEO belt (reference)', ...
    'Burn \DeltaV direction (length \propto |\DeltaV|)'}, ...
    'TextColor',[0.9 0.9 0.9], 'Box','off', 'FontSize',14, 'Location','northeast');

annotation('textbox', [0 0.92 1 0.06], 'String', ...
    'GEO End-of-Mission Disposal: Hohmann Transfer to Graveyard Orbit', ...
    'Color','w', 'FontSize',20, 'FontWeight','bold', 'EdgeColor','none', 'FitBoxToText','off', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
annotation('textbox', [0 0.875 1 0.05], 'String', ...
    sprintf('Total \\DeltaV = %.2f m/s     |     Transfer time = %.2f h     |     i = %g^\\circ', ...
            dVtot*1e3, TOF/3600, inc), ...
    'Color',[0.65 0.72 0.85], 'FontSize',14, 'EdgeColor','none', 'FitBoxToText','off', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
annotation('textbox', [0 0 1 0.04], 'String', ...
    sprintf(['Apsis radii exaggerated %g\\times from GEO for visibility (orbits drawn as exact ellipses ' ...
             'through those apsides); satellite not to scale'], k), ...
    'Color',[0.45 0.45 0.50], 'FontSize',12, 'FontAngle','italic', 'EdgeColor','none', ...
    'FitBoxToText','off', 'HorizontalAlignment','center', 'VerticalAlignment','middle');

% Save a high-resolution copy for slides (uncomment to use)
% exportgraphics(fig, 'geo_disposal_hohmann.png', 'Resolution',300, 'BackgroundColor','current');

%% Helper functions (MATLAB requires these at the end of a script)
function u = unitVec(v)
% Unit vector in the direction of v
u = v/norm(v);
end

function drawCone(tip, d, len, c)
% Solid 3D arrowhead: a cone with its point at "tip", pointing along d
d = unitVec(d(:));
E = null(d');                                    % two unit vectors perpendicular to d
[xc,yc,zc] = cylinder([0.4*len 0], 24);          % cone: base radius 0.4*len -> point
X = tip(1) + E(1,1)*xc + E(1,2)*yc + d(1)*len*(zc - 1);
Y = tip(2) + E(2,1)*xc + E(2,2)*yc + d(2)*len*(zc - 1);
Z = tip(3) + E(3,1)*xc + E(3,2)*yc + d(3)*len*(zc - 1);
surf(X, Y, Z, 'FaceColor',c, 'EdgeColor','none', 'AmbientStrength',0.7);
patch(X(1,:), Y(1,:), Z(1,:), c, 'EdgeColor','none', 'AmbientStrength',0.7);   % close the base
end