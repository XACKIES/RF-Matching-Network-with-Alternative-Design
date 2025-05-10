% ENE322 TransmissionLine : RF Matching Network Design project
% 65070502406 Kittiphop Phanthachart
% 65070502420 Than Thanyanothai
% 65070502498 Setthawut Kaweesukkaworakul



clear all;
clc;
ZL = 40+ j*15; % input impedance (Ohm)
Z0 = 30; % characteristic impedance (Ohm), must be real number
f0 = 8e9;
lambda0 = 3e8/f0;

% program starts here %
RL = real(ZL); % Resistance of the load
XL = imag(ZL); % Reactance of the load
Y0 = 1/Z0;


if ( RL > Z0 )
fprintf( '\n  Case ( RL > Z0 ) \n' );
B = (XL + (-1).^[1 0] * sqrt( RL/Z0)* sqrt( RL^2 + XL^2 -(Z0*RL)))/(RL^2 + XL^2)
X = (1./B) + (XL*Z0/RL) - (Z0./(B*RL))


elseif ( RL < Z0 )
fprintf('\n  Case ( RL < Z0 ) \n' );
B = (-1).^[1 0]* (sqrt( (Z0-RL)/RL))/Z0
X = (B.*Z0*RL)-XL

else 
fprintf('\n  Already Matching \n' );
end




nsol = length( B ); 
fprintf(1, '[L - Matching - Network] %d solution(s):', nsol );

for k=1:nsol

fprintf(1, '\nSolution #%d\n', k );

if (B(k)>0 && X(k)>0)
fprintf(1, '\n  Series L Parallel C   \n');
L(k) = X(k)/(2*pi*f0);
C(k) = B(k)/(2*pi*f0);

elseif (B(k)<0 && X(k)<0)
fprintf(1, '\n  Series C Parallel L   \n');
L(k) = -1/(2*pi*B(k)*f0);
C(k) = -1/(2*pi*X(k)*f0);
end

fprintf(1, ' ----[B = %g]---- \n', B(k) );
fprintf(1, ' ----[X = %g]---- \n', X(k) );
fprintf(1, ' L = %g\n', L(k) );
fprintf(1, ' C = %g\n', C(k) );

end



f = linspace(0, 15e9, 1000); 

lambda = 3e8 ./ f; 

Gamma = zeros(nsol, length(f));

for k = 1:nsol

    L_load = XL/(2*pi.*f0);
    Za     = RL + j*2*pi.*f*L_load;            % series L

    if B(k)>0 && X(k)>0
        Zb      = 1./( j*2*pi.*f * C(k) );     % shunt C
        Zc      = ( Za.^-1 + Zb.^-1 ).^-1;      % parallel
        Z_total = Zc + j*2*pi.*f*L(k);         % series L


    elseif B(k)<0 && X(k)<0
        Zb      = j*2*pi.*f*L(k);              % shunt L
        Zc      = ( Za.^-1 + Zb.^-1 ).^-1;      % parallel
        Z_total = Zc + 1./( j*2*pi.*f * C(k) ); % series C
    end

    Gamma(k,:) = abs( (Z_total - Z0) ./ (Z_total + Z0) );
end

figure;
plot(f / 1e9, Gamma(1,:), 'LineWidth', 2); hold on;
plot(f / 1e9, Gamma(2,:), 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (|\Gamma|)');
title('L- Network Tunner ');
legend(' Series C Parallel L   Solution 1', ' Series L Parallel C  Solution 2');
grid on;

target = 0.2;
for k = 1:size(Gamma,1)
    G = Gamma(k,:);
    % down‐crossing and up‐crossing
    idx_down = find( G(1:end-1)> target  & G(2:end)<= target, 1, 'first');
    idx_up   = find( G(1:end-1)<=target  & G(2:end)>  target, 1, 'last');

    if isempty(idx_down)
      fL = 0;
    else
      fL = interp1( G(idx_down:idx_down+1), f(idx_down:idx_down+1), target );
    end

   if isempty(idx_up)
      fH = 0;
    else
      fH = interp1( G(idx_up:idx_up+1), f(idx_up:idx_up+1), target );
    end

    BW  = fH - fL;
    FBW = BW / f0;

    fprintf('Solution %d:\n', k);
    fprintf('  f_L = %.2f GHz\n', fL/1e9);
    fprintf('  f_H = %.2f GHz\n', fH/1e9);
    fprintf('  BW  = %.2f GHz\n', BW/1e9);
    fprintf('  FBW = %.2f %%\n\n', FBW*100);
end
