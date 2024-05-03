function NMSE = THUNDERR_NMSE(X,ti,dt,tf,xp,yp,vrt,alphart,fs,kv,kalpha)

%% ----------------------Thunderr_02_NMSE____------------------------------
% Thunderr_02_NMSE function is an advanced analitical model for the
% reconstruction of the horizontal mean wind velocity that accours during a
% downburst (thunderstorm)... contina, spiega cosa è NMSE


%% -----------------Model external parameters-----------------------------%
xc0 = X(1);        % - x  component of  downburst touch down (m)
yc0 = X(2);        % - y  component of  downburst touch down (m)
R = X(3);          % - downburst downdraft radius (m)
psi = X(4);        % - psi = Rmax/R Radius of maximum radial velocity 
v1 =X(5);          % - downburst maximum radial velocity (m/s), Holmes and Oliver
                   %   function (2000)
T1 = X(6);         % - period of linear intensification of the downburst (s), Chay at al. function (2006)
Tf = X(7);         % - total duration of the downburst (s), Chay at al. function (2006)
v2 = X(8);         % - downburst translation velocity (m/s)
beta2 = X(9);      % - downburst translation direction (polar angle from the East) (deg)
v3 = X(10);        % - ABL background wind speed (m/s)
beta3 = X(11);     % - ABL background wind direction (polar angle from the East) (deg)
%-------------------------------------------------------------------------%

%% -------------------Model Necessary Parameters -------------------------%
% ti       - initial simulation time (s), usually = 0 s
% dt       - simulation time step (s)
% tf       - final simulation time (s)
% xp       - x station location component 
% yp       - y station location component
% vxrt     - moving average(T_ma = 30 s) recorded (r) data along x - axis
% vyrt     - moving average(T_ma = 30 s) recorded (r) data along y - axis
% vrt      - moving average(T_ma = 30 s) wind speed recorded (r) data 
% alphart  - moving average(T_ma = 30 s) wind direction recorded (r) data
% fs       - recorded data sampling frequency (10 Hz mostly)

% kv       - wheight cofficient NMSE_vy,    0 <= kv <= 1
% kalpha   - wheight cofficient NMSE_alpha, 0 <= kalpha <= 1
%-------------------------------------------------------------------------%

%% -------------------Model Internal Parameters----------------------------
% Internal parameters -  Holmes and Oliver function 2000
Rs = R;                  % Exponential decay constant (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rmax = psi*R;             % Radius of maximum radial velocity (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Internal parameters -  Auxiliary function translational component
RE = 2*Rmax;             % radius of maximum outflow extension (m)
a = 0.2;                      % smooth factor
RE_max = (1 + a)*RE;
%-------------------------------------------------------------------------%
% Internal parameters -  Intensity decay function translational component
b = 0.2;                 % smooth factor         
Ti = b.*T1;                    
Tf_max = (1 + b).*Tf;  
%-------------------------------------------------------------------------%
%% ------------------ Starting of Thunderr_02_NMSE------------------------%
% Time axis
t = (ti:dt:tf-dt)';
nt = dt*fs;

%% ----------------------(1) - Stationary downburst-----------------------%

% Holmes and Oliver function 2000

vr = @(s) (v1.*s./Rmax).*(s <= Rmax)...
                      + (v1.*exp(-(((s-Rmax)./Rs).^2))).*(s > Rmax);
% Note
% s is a dummy variable
%-------------------------------------------------------------------------%

% Chay at al. 2006 - intensity decat function radial component
% Enssuring that Tf is grater than T1
if Tf < T1           
    error (' Tf <= T1, Tf shoud be > T1')
end
% Exponential decay time constant (sec)
c = (Tf - T1)./(log(10));         
PIr = @(s) (s./T1).*(s<=T1) + (exp(-((s-T1)./c))).*(s>T1);
%-------------------------------------------------------------------------%

%% ----------------------(2) - Translating downburst----------------------%

% Downburst path
xct = xc0 + v2*cosd(beta2)*t;
yct = yc0 + v2*sind(beta2)*t;
%-------------------------------------------------------------------------%
% Relative distance between the downburst center and the station
rt = sqrt((xp - xct).^2 + (yp - yct).^2);

% Relative direction between the downburst center and the station
% (polar angle i.e. from East)
beta1t = atan2d((yp - yct),(xp - xct));
%-------------------------------------------------------------------------%
% Xhelaj et al. (2020) Auxiliary function translational component

Delta = @(s) (1).*((0 <= s)&(s <= RE))...
              + ((1/2).*(1 + cos((pi/(a*RE)).*(s - RE)))).*((RE < s)&(s <= RE_max))...
              + (0).*(s > RE_max);
%-------------------------------------------------------------------------%
% Xhelaj et al. (2020) Intensity - decay function translational component

PIt =  @(s) ((1/2).*(1 + cos((pi./Ti).*(s - Ti)))).*((0 <= s)&(s <= Ti))...
          + (1).*((Ti < s)&(s <= Tf))...
          + ((1/2).*(1 + cos((pi./(b.*Tf)).*(s - Tf)))).*((Tf < s)&(s <= Tf_max))...
          + (0).*(s > Tf_max);

%% ------------------------ Simulated data---------------------------------
vxt = vr(rt).*PIr(t).*cosd(beta1t) + v2.*Delta(rt).*PIt(t).*cosd(beta2) + v3*cosd(beta3);
vyt = vr(rt).*PIr(t).*sind(beta1t) + v2.*Delta(rt).*PIt(t).*sind(beta2) + v3*sind(beta3);

vt = sqrt(vxt.^2 + vyt.^2);
betat = atan2d(vyt,vxt);
betat = mod(betat,360);
alphat = mod(270 - betat,360);
%--------------------------------------------------------------------------
%% -------------------------- Recorded data--------------------------------

vrt = vrt(1:nt:end);
alphart = alphart(1:nt:end);
%--------------------------------------------------------------------------
%% --------------------- Shifting the  simulated data ---------------------
[~,ind1] = max(vrt);
[~,ind2] = max(vt);
Shift_abs = abs(ind2 - ind1);
if  ind1 > ind2
    Shift = Shift_abs;
else
    Shift = - Shift_abs;
end

vt = circshift(vt,Shift);
alphat = circshift(alphat,Shift);
%--------------------------------------------------------------------------
%% ------------------ Objective function evaluation------------------------

p = 2;                                      % generalized vector p-norm.

NMSE_v = norm(vt - vrt,p)/(norm(vrt,p));
% NOTE: CIRCULAR STATISTICS HERE!

Xt = [sind(alphat), cosd(alphat)];
Xrt = [sind(alphart), cosd(alphart)];
NMSE_alpha = norm(Xt - Xrt,'fro')/(norm(Xrt,'fro'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NMSE =  kv*NMSE_v +kalpha*NMSE_alpha;                 % OBJECTIVE FUNCTION  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end