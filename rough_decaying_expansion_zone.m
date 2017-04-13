% Evalute a possible curve to find the extent of the shear zone
% EJR 2017
%
%  Use the assumption that a horizontal strip of finite L0 extent becomes 
% a cosine with a range 0 to pi, and matches plate dispacement on one 
% side and zero movement on the other. 
%  Let this strip describe the shear zone near the plate edge.
%  Assume strips higher up have the same form of dispacement curve, 
% but a larges horizontal extent to accommodate dilatancy of sublayer(s).
%
% NOTES
% 1. Check: almost surely need to put (2/zp) back into elongation equation.
%
% 2. Dilatancy is modelled as eV_actual = eM*(1-exp(-gamma/gamma0))
%
% 3. The internal energy of each strip is then calculated. The model
%   assumes internal energy = gravitational energy due to uplift, plus
%   'p dV' work due to dilatancy, plus a term for elastic shear 'strain'
%   work which seems to be needed for plausible shear zones to be 
%   obtained, but which needs a stronger physical justification if possible
%
% 4. The draft version starts with 'small zp' so a strip depth ~ uniform 
%    For evaluating forces with variable zp, we MIGHT need to account for 
%   the depth of many parts of the shear zone being smaller than implied?

% PART 1: Evaluate one perimeter curve to find extent of shear zone

L0 = 5;      % mm, horizontal extent of shear zone at base strip
zStep = 1;   % mm, Strip thickness (should be 'small' c.f. zp)
nSteps = 100; % number of strips to top of modelled region

H  = nSteps*zStep;  % mm, Depth of ballotini to initial plate position
eM = 0.05;   % maximal dilatancy (volume strain) of material

zp = 4;      % mm, plate displacement
gam0 = 0.12; % radians, scale factor for material dilatancy

rho = 1700;  % kg.m^3, bulk density - irrelevant, cancels out
K   = 1;     % Arbitary number to produce 'effective elastic modulus'

listL = zeros(nSteps,1); % list to store strip widths
listU = zeros(nSteps,1); % list to store strip energies
L = L0;      % initialise strip width for Euler ODE method
for lp = 1:nSteps
	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ) * (zStep*2/zp);
	
	listL(lp) = L;
	L = L+k*zStep;
	
	Zn = H - (lp-1)*zStep; % depth of strip base
	enGPE = zp*L/2;        % estimate gravitational internal energy added
	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
	
	listU(lp) = enGPE + enPDV + enSHR;
end

listAltitudes = 0:zStep:(nSteps-1)*zStep;
figure(1)
plot(listL, listAltitudes)
sum(listU)
xlim([0 100])
ylim([0 100])
xlabel('x-position from plate edge/ mm')
ylabel('height above plate start-position / mm')
title('perimeter of shear zone for a given L0')
set(gca, 'fontSize', 12)

%% PART 2
%  Evaluate a family of curves to find the extent of the shear zone
%  Consider several possible values of L0. 
%  Choose the one that produces the shear zone with the lowest energy
% 
% NOTES:
% 1. Run PART 1 first to initialise system parameters

nLs = 20;
listsL = zeros(nSteps, nLs); % array to store several perimeter curves
listsUT = zeros(nLs,1);      % array to store total energy of each zone
close(1) % Cleanly overwrite plot from PART 1
listL0s = 1:1*1:(nLs*1);

for lp2 = 1:nLs;
L0 = listL0s(lp2)  % mm
zStep = 1;         % mm
% nSteps = 30; 

eM = 0.05; 

zp = 4; % mm
gam0 = 0.12; % radians

listL = zeros(nSteps,1);
L = L0;

for lp = 1:nSteps
	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ) * (zStep*2/zp);
	
	listL(lp) = L;
	L = L+k*zStep;
	
	Zn = H - (lp-1)*zStep; % depth of strip base
	enGPE = zp*L/2;        % estimate gravitational internal energy added
	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
	
	listU(lp) = enGPE + enPDV + enSHR;
end
% for lp = 1:nSteps
% 	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ); %  * (2/zp);
% 	
% 	listL(lp) = L;
% 	L = L+k*zStep;
% 	
% 	
% 	Zn = H - (lp-1); % depth of strip base
% 	enGPE = zp*L/2;
% 	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
% 	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
% 	
% 	listU(lp) = enGPE + enPDV + enSHR;
% end

listsL(:,lp2) = listL;
listsUT(lp2)  = sum(listU);
figure(1)
hold on
plot(listL, 1:nSteps)
end

figure(1)
hold off
title('Perimeter of possible shear zones for various base extents, L0')
xlabel('x-position from plate edge/ mm')
ylabel('height above plate start-position / mm')
set(gca, 'fontSize', 12)
xlim([0 100])
ylim([0 100])

figure(2)
plot(listL0s, listsUT);
title('total energy of shear zone given base extent L0')
xlabel('Horizontal extent of shear zone base, L0 / mm')
ylabel('Total energy of shear zone / arb')
set(gca, 'fontSize', 12)

%% PART 3
%  Now loop through successive z_plate values
%  Evaluate the (minimum possible, hence equilibrium) shear zone energy
% for various plate alititudes. 
%  Hence equate shear zone energy (and uplift column energy) to 
% identify a corresponding force. 
%  Note that the real peak force may be limited by void collapse which 
% isn't included in this model - could be introduced as a cut-off. 

nZs = 20;
Uz = zeros(nZs,1);
figure(1)

for lp3 = 1:nZs
zp = lp3;

nLs = 20;
listsL = zeros(nSteps, nLs);
listsUT = zeros(nLs,1);

% LOOP -> loop
for lp2 = 1:nLs;
L0 = listL0s(lp2)  % mm
zStep = 1;         % mm
% nSteps = 30; 

eM = 0.05; 

% zp = 4; % mm
gam0 = 0.12; % radians

listL = zeros(nSteps,1);
L = L0;

for lp = 1:nSteps
	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ) * (zStep*2/zp);
	
	listL(lp) = L;
	L = L+k*zStep;
	
	Zn = H - (lp-1)*zStep; % depth of strip base
	enGPE = zp*L/2;        % estimate gravitational internal energy added
	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
	
	listU(lp) = enGPE + enPDV + enSHR;
end
% for lp = 1:nSteps
% 	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ); %  * (2/zp);
% 	
% 	listL(lp) = L;
% 	L = L+k*zStep;
% 	
% 	
% 	Zn = H - (lp-1); % depth of strip base
% 	enGPE = zp*L/2;
% 	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
% 	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
% 	
% 	listU(lp) = enGPE + enPDV + enSHR;
% end

listsL(:,lp2) = listL;
listsUT(lp2)  = sum(listU);
% figure(1)
% hold on
% plot(listL, 1:nSteps)
end

Uz(lp3) = min(listsUT);
end

figure(3)
plot(Uz)
title('sand internal energy')
xlabel('z plate ')
ylabel('internal energy / very arb units')

% calculate column weight in same arb units to estimate breakout factor
% GPE of one strip (H strips) = zp * L/2 in arb units. Both L, zp in mm. 
% wColumn = W_real / (rho g L delH?)
% wColumn = H*(D/2)
% wColumn*zp = H*D/2*zp
% say column half width = 20 mm, H = 100 mm
wColumn = 100*40;

delU = Uz(2:end) - Uz(1:(end-1));
figure(4)
plot(delU + wColumn)
title('energy change with plate step')
xlabel('z plate ')
ylabel('force (very roughly) / arb')
% ylim([0 7000])


%%
nZs = 20;
Uz = zeros(nZs,1);
figure(1)

for lp3 = 1:nZs
zp = lp3;

nLs = 20;
listsL = zeros(nSteps, nLs);
listsUT = zeros(nLs,1);
% close(1)

for lp2 = 1:nLs;
L0 = lp2*2; % mm
zStep = 1; % mm
nSteps = 30;

eM = 0.05; 

% zp = 5; % mm
gam0 = 0.12; % radians

listL = zeros(nSteps,1);
L = L0;
for lp = 1:nSteps
	k = L* eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0)) ); %  * (2/zp);
	
	listL(lp) = L;
	L = L+k*zStep;
	
	
	Zn = max( H - (lp-1) - zp ,  0); % depth of strip base. ASSUME TOP CLEARs
	enGPE = zp*L/2;
	enPDV = Zn*L*eM*(1 - (6.209/pi)*exp(-(zp*pi)/(2*L*gam0) ) );
	enSHR = (1/16)*K*Zn*pi^2*zp^2/L;
	
	listU(lp) = enGPE + enPDV + enSHR;
end

listsL(:,lp2) = listL;
listsUT(lp2)  = sum(listU);
figure(1)
hold on
plot(listL, 1:nSteps)
end

Uz(lp3) = min(listsUT);

end

figure(3)
plot(Uz)
title('sand internal energy')
xlabel('z plate ')
ylabel('internal energy / very arb units')

% calculate column weight in same arb units to estimate breakout factor
% GPE of one strip (H strips) = zp * L/2 in arb units. Both L, zp in mm. 
% wColumn = W_real / (rho g L delH?)
% wColumn = H*(D/2)
% wColumn*zp = H*D/2*zp
% say column half width = 20 mm, H = 100 mm
wColumn = 100*40;

delU = Uz(2:end) - Uz(1:(end-1));
figure(4)
plot(delU)
title('energy change with plate step')
xlabel('z plate ')
ylabel('force (very roughly) / arb')
% ylim([0 7000])
