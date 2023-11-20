clear all;
close all;

Fe = 6000;          % Fréquence d'échantillonage
Rb = 3000;          % Débit binaire de la transmission
N_bits = 300000;    % Nombre de bits à transmettre
f_p = 2000;         % Fréquence porteuse
M = 8;              % Ordre de la modulation
EbsurN0 = 0:1:6;    % Rapports signal sur bruit à tester
alpha = 0.20;       % Roll-off
L = 10;             % Span

Te = 1 / Fe;                        % Période d'échantillonage
Ts = 1 / Rb * log2(M);              % Période par symbole
Ns = fix(Ts / Te);                  % Facteur de suréchantillonage
bits = randi([0, 1], 1, N_bits);    % Bits à transmettre
t = (0:Ns*N_bits/2-1) * (1/Fe);     % Échelle temporelle

mat_kron = [1 zeros(1, Ns-1)];
h = rcosdesign(alpha, L, Ns);
Nt = L / 2;
retard = Nt * Ns;

%% Génération du signal
bits_triplets = reshape(bits, 3, N_bits / 3);                           % Mapping
symboles = bi2de(bits_triplets.');                                      %
signal = pskmod(symboles, 8, 0, 'gray', 'PlotConstellation', true).';   %
signal_surech = kron(signal, mat_kron);                                 % Suréchantillonnage
signal_surech_retard = cat(2, signal_surech, zeros(1, retard));         % Gestion du retard
x_e = filter(h, 1, signal_surech_retard);                               % Filtrage
x_e = x_e(1, retard+1 : end);                                           % Gestion du retard

teb = zeros(1, length(EbsurN0));
for i=1:length(EbsurN0)
    %% Canal de propagation AWGN
    M = 8;                                                  % Ordre de la modulation
    Px_e = mean(abs(x_e).^2);                               % Puissance du signal transmis
    sigmacarre = Px_e*Ns / (2*log2(M)*10^(EbsurN0(i)/10));  % Puissance du bruit pour le rapport Eb/N0 souhaité
    bruit_I = sqrt(sigmacarre) * randn(1, length(x_e));     % Bruit réel
    bruit_Q = sqrt(sigmacarre) * randn(1, length(x_e));     % Bruit complexe
    bruit = bruit_I + 1i*bruit_Q;                           % Bruit
    y = x_e + bruit;                                        % Signal bruité

    %% Filtre de réception
    y_protection_retard = cat(2, y, zeros(1, retard));  % Gestion du retard
    z = filter(h, 1, y_protection_retard);              % Retour en bande de base
    z = z(1, retard+1 : end);                           % Gestion du retard

    %% Échantillonnage
    n0 = 1;                             % Instant d'échantillonnage (retard géré avant)
    z_echantillonne = z(n0 : Ns : end); % Échantillonnage à n0+m*Ns
    
    %% Plot constellation
    scatterplot(z_echantillonne);
    xlabel("En phase");
    ylabel("En quadrature");
    title("Constellation en sortie de l'échantillonneur (Eb/N0=" + EbsurN0(i) + "dB)");

    %% Décision et démapping
    demod_symboles = pskdemod(z_echantillonne, 8, 0, 'gray')';
    demod_bits = de2bi(demod_symboles)';
    demod_bits = reshape(demod_bits, 1, N_bits);

    %% Taux d'erreur binaire
    teb(i) = 1-sum(demod_bits == bits) / N_bits;
    
end

DSP_x_e = pwelch(x_e, [], [], [], Fe, 'twosided');                      % DSP estimée de x_e   
f = linspace(-Fe/2, Fe/2, length (DSP_x_e));                            % Échelle fréquentielle
teb_th = 2*qfunc(sqrt(2*log2(M)*10.^(EbsurN0/10))*sin(pi/M))/log2(M);   % TEB théorique pour une chaine M-PSK

%% Figures
figure;
    nexttile;   % Tracé TEB
    semilogy(EbsurN0, teb);
    xlabel("Eb/N0 (dB)");
    ylabel("TEB");
    title("TEB en fonction du rapport signal à bruit par bit");

figure;
    nexttile;   % Tracé TEB simulé vs théorique
    semilogy (EbsurN0, teb);
    hold on;
    semilogy (EbsurN0, teb_th);
    hold off;
    xlabel("Eb/N0 (dB)");
    ylabel("TEB");
    legend ("TEB simulé", "TEB théorique");
    title("Comparaison des TEB simulé et théorique");

figure;
    nexttile;   % Tracé DSP x_e
    semilogy(f, fftshift(DSP_x_e));
    xlabel("Fréquence (Hz)");
    ylabel("DSP");
    title("DSP de l'enveloppe convexe associée au signal transmis sur fréquence porteuse");
