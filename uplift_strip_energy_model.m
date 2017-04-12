zp = 0.012;  % plate displacement
epsM = 0.05; % Max volume or area dilation
H = 0.2;     % ballotini depth / m
gam = 0.14;   % gamma_zero for dilatancy correlation
mu = 0.6;    % Co-efficient of friction to bodge a shear modulus

xs = 1E-4:1E-4:0.1; % Horizontal extent of shear strip.

Up_GPE = xs.*(zp/2);
Up_Exp = max( H*epsM.*xs.*( 1 - (6.21/pi).*exp(-(pi*zp)./(2.*xs*gam)) ),0);
Up_Shr = mu.*(1./xs).*(H*pi^2*zp^2)/16;

Up = Up_GPE + Up_Exp + Up_Shr;

figure(2)
plot(xs, Up, 'k', 'linewidth',2)
hold on
 plot(xs, Up_Exp, 'b', 'linewidth',2)
 plot(xs, Up_Shr, 'r', 'linewidth',2)
 plot(xs, Up_GPE, 'g', 'linewidth',2)
hold off
legend('total energy', 'expansion energy', 'shear work', 'gravity')
ylim([0 0.002])
ylabel('Energy / arb units')
xlabel('Extent of shear zone / m')
set(gca, 'fontsize', 12)

%%
XS = 0.015
xx0 = 0:1E-4:0.1;
yy0 = zeros(length(xx0), 1);
yy1 = yy0;
yy1(xx0<=XS) = (zp/2)*(1+cos((pi*xx0(xx0<=XS) /XS)));

figure(2)
plot(xx0,yy1)


%% Try estimating shear zone width, allowing for expansion of lower layers
% by including them as an increase in the effective plate movement to be 
% accommodated
zp = 0.004;
figure(3)
hold on
z1 = 0.25;
z2 = 0.05;
zstep = 0.05;
listH = [z1:-zstep:z2];
for lp = 1:5;

H =  listH(lp)
if(lp>1)
	lastXS = xs(Up==min(Up));
	gamMax = (1/lastXS)* 1 * (zp/2);
	espV = epsM*(1 - exp(-(gamMax/gam)) )
	zp = zp+(zstep)*espV*0.5
end

Up_GPE = xs.*(zp/2);
Up_Exp = max( H*epsM.*xs.*( 1 - (6.21/pi).*exp(-(pi*zp)./(2.*xs*gam)) ),0);
Up_Shr = mu.*(1./xs).*(H*pi^2*zp^2)/16;

Up = Up_GPE + Up_Exp + Up_Shr;

plot(xs,Up)
end
% plot(xs, Up,'r')
legend('1','2','3','4','5','6')
H
hold off
ylim([0 0.004])
