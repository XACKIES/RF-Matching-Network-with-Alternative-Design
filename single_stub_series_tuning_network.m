% ENE322 TransmissionLine : RF Matching Network Design project
% 65070502406 Kittiphop Phanthachart
% 65070502420 Than Thanyanothai
% 65070502498 Setthawut Kaweesukkaworakul


clear all;
clc;
ZL = 40 + j*15; % input impedance (Ohm)
Z0 = 30; % characteristic impedance (Ohm), must be real number
f0 = 8e9;
lambda0 = 3e8/f0;

% program starts here %
YL = 1/ZL; % load ad
GL = real(YL); % conductance of the load
BL = imag(YL); % susceptance of the load
Y0 = 1/Z0; % characteristic admittance of the line

% obtain t
if ( GL ~= Y0 )
% two solution. Put them in a vector 't'
t = ( BL + (-1).^[1 0] * sqrt( GL*( (Y0-GL)^2 + BL^2 )/Y0 ) ) / ( GL - Y0 );
else
% one solution
t = -BL/(2*Y0);
end

% obtain X
X = ( GL^2*t - (Y0-BL*t).*(BL + Y0*t) ) ./ ( Y0*(GL^2 + (BL + Y0*t).^2 ) );

% obtain normalized length, norm_ls = ls/lambda, of the short-circuit stub

norm_ls = -atan( X/Z0 ) / (2*pi); 
norm_ls( norm_ls < 0 ) = norm_ls( norm_ls < 0 ) + 1/2;
ls = norm_ls * lambda0;

% obtain normalized length, norm_lo = lo/lambda, of the open-circuit stub
norm_lo = atan( Z0./X ) / (2*pi);
norm_lo( norm_lo < 0 ) = norm_lo( norm_lo < 0 ) + 1/2;
lo = norm_lo * lambda0;

% obtain the normalized distance, norm_d = d/lambda, of the stub
% Note that t can be a vector, so we can have multiple solutions of norm_d
norm_d = atan( t ) / (2*pi);
norm_d( t<0 ) = norm_d( t< 0 ) + 1/2;
d = norm_d * lambda0;

% Print out the solutionsf
nsol = length( norm_d ); % number of solution
fprintf(1, '[Single-stub series tuner] %d solution(s):', nsol );
for k=1:nsol
fprintf(1, '\nSolution #%d\n', k );
fprintf(1, ' Distance of the stub: d/lambda = %g\n', norm_d(k) );
fprintf(1, ' d= %g\n', d(k) );
fprintf(1, ' Short circuit: ls/lambda = %g\n', norm_ls(k) );
fprintf(1, ' ls = %g\n', ls(k) );
fprintf(1, ' Open circuit: lo/lambda = %g\n', norm_lo(k) );
fprintf(1, ' lo = %g\n', lo(k) );

end


f = linspace(0, 12e9, 1000); 
lambda = 3e8 ./ f;
beta = 2 * pi ./ lambda;
Gamma_short = zeros(nsol, length(f));


for i = 1:nsol
    

    Zin = Z0 .* (ZL + 1j * Z0 .* tan(beta .* d(i))) ./ (Z0 + 1j * ZL .* tan(beta .* d(i)));
    Z_stub = (1)*j * Z0 .* tan(beta * ls(i));

    Z_total = Zin + Z_stub;
    Gamma_short(i, :) = abs((Z_total - Z0) ./ (Z_total + Z0));
end

figure;
plot(f / 1e9, Gamma_short(1,:), 'LineWidth', 2); hold on;
plot(f / 1e9, Gamma_short(2,:), 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (|\Gamma|)');
title('Single Stub series tunning (Short)');
legend('Short Solution 1', 'Short Solution 2');

grid on;

%-----------------------------------------------------------------------------------------------------
target = 0.2;
for k = 1:size(Gamma_short,1)
    G = Gamma_short(k,:);
    
    idxL = find(  f(1:end-1)<f0  &  G(1:end-1)>target  &  G(2:end)<= target, 1, 'last' );
    if isempty(idxL)
      fL = NaN;
    else
      fL = interp1( G(idxL:idxL+1), f(idxL:idxL+1), target );
    end

    
    idxH = find(  f(2:end)>f0    &  G(1:end-1)<=target  &  G(2:end)> target,  1, 'first' );
    if isempty(idxH)
      fH = NaN;
    else
     
      fH = interp1( G(idxH:idxH+1), f(idxH:idxH+1), target );
    end

    BW  = fH - fL;
    FBW = BW/f0*100;

    fprintf('\n Short Solution %d near f0:\n', k);
    fprintf('  f_L = %.3f GHz\n', fL/1e9);
    fprintf('  f_H = %.3f GHz\n', fH/1e9);
    fprintf('  BW  = %.3f GHz\n', BW/1e9);
    fprintf('  FBW = %.2f %%\n\n', FBW);
end





Gamma_open = zeros(nsol, length(f));
for i = 1:nsol

    Zin = Z0 .* (ZL + 1j * Z0 .* tan(beta .* d(i))) ./ (Z0 + 1j * ZL .* tan(beta .* d(i)));
    Z_stub = (-1)*j * Z0 .* cot(beta * lo(i));

    Z_total = Zin + Z_stub;
    Gamma_open(i, :) = abs((Z_total - Z0) ./ (Z_total + Z0));
end

figure;
plot(f / 1e9, Gamma_open(1,:), 'LineWidth', 2); hold on;
plot(f / 1e9, Gamma_open(2,:), 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (|\Gamma|)');
title('Single Stub series tunning (Open) ');
legend('Open Solution 1', 'Open Solution 2');

grid on;


target = 0.2;
for k = 1:size(Gamma_open,1)
    G = Gamma_open(k,:);
    
    idxL = find(  f(1:end-1)<f0  &  G(1:end-1)>target  &  G(2:end)<= target,  1, 'last' );
    if isempty(idxL)
      fL = NaN;
    else
      fL = interp1( G(idxL:idxL+1), f(idxL:idxL+1), target );
    end


    idxH = find(  f(2:end)>f0    &  G(1:end-1)<=target  &  G(2:end)> target,  1, 'first' );
    if isempty(idxH)
      fH = NaN;
    else

      fH = interp1( G(idxH:idxH+1), f(idxH:idxH+1), target );
    end

    BW  = fH - fL;
    FBW = BW/f0*100;

    fprintf('\n Open Solution %d near f0:\n', k);
    fprintf('  f_L = %.3f GHz\n', fL/1e9);
    fprintf('  f_H = %.3f GHz\n', fH/1e9);
    fprintf('  BW  = %.3f GHz\n', BW/1e9);
    fprintf('  FBW = %.2f %%\n\n', FBW);
end
