% A plotting script

% dfigure.mat, initializes a docked figure. cf. utils repository.
% Lineout_Radial.m, makes a radial lineout, cf. procVMI repository.

%% plot Abel inversion results

PtoE = 0.0887; % conversion from pixel to energy, [eV/pixel^2]

dfigure;
imagesc(AngIntegrated); title('cross section');

dfigure;
imagesc(Ring); title('3D density');

% lineout
ltInv = Lineout_Radial(AngIntegrated,Centre,12); % Centre = [row center, coloum center].
% plot(PtoE*(1:x).^2, squeeze(lnt(3,:)+lnt(8,:)));
x = size(ltInv,2);
dfigure;
% Shift by half a pixel for pixel-centered data. 
% Data is function of pixel, so dividing by 1:x to get density. 
plot(PtoE*((1:x)-0.5).^2, ltInv(3,1:x)./((1:x)-0.5));


