clear all;
close all;

Fe = 24000;         % Fréquence d'échantillonage
Rb = 3000;          % Débit binaire
N_bits = 10000;     % Nombre de bits à transmettre
f_p = 2000;         % Fréquence porteuse
M = 4;              % Ordre de la modulation
EbsurN0 = 0:1:6;    % Rapports signal sur bruit à tester
alpha = 0.35;       % Roll off
L = 10;             % Span

Te = 1 / Fe;                        % Période d'échantillonage
Ts = 1 / Rb * log2(M);              % Période par symbole
Ns = fix(Ts / Te);                  % Facteur de suréchantillonage
bits = randi([0, 1], 1, N_bits);    % Bits à transmettre
t = (0:Ns*N_bits/2-1) * (1/Fe);     % Échelle temporelle

%% Séparation bits pairs et impairs
bits_pairs = zeros(1, N_bits / 2);
bits_impairs = zeros(1, N_bits / 2);
for i = 1:(N_bits/2)
    bits_pairs(i) = bits(2*i);
    bits_impairs(i) = bits(2*i-1);
end

%% Génération du signal
mat_kron = [1 zeros(1, Ns-1)];
h = rcosdesign(alpha, L, Ns);
Nt = L / 2;
retard = Nt * Ns;

%% ... sur la voie en phase (bits pairs)
symboles_pairs = (bits_pairs==1) - (bits_pairs==0);                     % Mapping
Mod_pairs = kron(symboles_pairs, mat_kron);                             % Suréchantillonage
Mod_pairs_protection_retard = cat(2, Mod_pairs, zeros(1, retard));      % Gestion du retard
x_I = filter(h, 1, Mod_pairs_protection_retard);                        % Filtrage
x_I = x_I(1, retard+1 : end);                                           % Gestion du retard

%% ... sur la voie en quadrature (bits impairs)
symboles_impairs = (bits_impairs==1) - (bits_impairs==0);               % Mapping
Mod_impairs = kron(symboles_impairs, mat_kron);                         % Suréchantillonage
Mod_impairs_protection_retard = cat(2, Mod_impairs, zeros(1, retard));  % Gestion du retard
x_Q = filter(h, 1, Mod_impairs_protection_retard);                      % Filtrage
x_Q = x_Q(1, retard+1 : end);                                           % Gestion du retard

%% Transposition de fréquence
x_e = x_I + 1i*x_Q;                     % Enveloppe convexe associée au signal
x = real(exp(2*1i*pi*f_p*t) .* x_e);    % Signal modulé sur porteuse

teb = zeros(1, length(EbsurN0));
for i=1:length(EbsurN0)
    %% Canal de propagation AWGN
    Px = mean(abs(x).^2);                                   % Puissance du signal transmis
    sigmacarre = Px*Ns / (2*log2(M)*10.^(EbsurN0(i)/10));   % Puissance du bruit pour le rapport Eb/N0 souhaité
    bruit = sqrt(sigmacarre) * randn(1, length(x));         % Bruit
    y = x + bruit;                                          % Signal bruité

    %% Filtre de réception
    x_I_estime = cos(2*pi*f_p*t) .* y;                                  % Reconstruction du signal x_I
    x_Q_estime = sin(2*pi*f_p*t) .* y;                                  % Reconstruction du signal x_Q
    x_estime = x_I_estime - 1i*x_Q_estime;                              % Reconstruction du signal x
    x_estime_protection_retard = cat(2, x_estime, zeros(1, retard));    % Gestion du retard
    z = filter(h, 1, x_estime_protection_retard);                       % Retour en bande de base
    z = z(1, retard+1 : end);                                           % Gestion du retard

    %% Échantillonnage
    n0 = 1;                                 % Instant d'échantillonnage (retard géré avant)
    z_echantillonne = z(n0 : Ns : end);     % Échantillonnage à n0+m*Ns
    z_I = real(z_echantillonne);            % Partie réelle
    z_Q = imag(z_echantillonne);            % Partie imaginaire

    %% Décision et démapping
    bits_pairs_demap = 1*(z_I >= 0) + 0*(z_I < 0);      % Démapping pour les bits pairs
    bits_impairs_demap = 1*(z_Q >= 0) + 0*(z_Q < 0);    % Démapping pour les bits impaires
    bits_demap = zeros(1, N_bits);                      % Bits
    for j=1:N_bits/2
        bits_demap(2*j) = bits_pairs_demap(j);
        bits_demap(2*j-1) = bits_impairs_demap(j);
    end

    %% Taux d'erreur binaire
    teb(i) = 1-sum(bits_demap == bits) / N_bits;
end

DSP_x_I = pwelch(x_I, [], [], [], Fe, 'twosided');  % DSP estimée de x_I
DSP_x_Q = pwelch(x_I, [], [], [], Fe, 'twosided');  % DSP estimée de x_Q
DSP_x = pwelch(x, [], [], [], Fe, 'twosided');      % DSP estimée de x
teb_th = 2*qfunc(sqrt(2*log2(M)*10.^(EbsurN0/10))*sin(pi/M))/log2(M);
f = linspace (-Fe/2, Fe/2, length (DSP_x));         % Échelle fréquentielle

%% Figures
% figure;
%     nexttile;   % Tracé diagramme de l'oeil
%     plot(reshape(z, Ns, length(z)/Ns));
figure;
    nexttile;   % Tracé signal x_I
    plot(t, x_I);
    xlim([0 t(length(t))]);
    xlabel("Temps (s)");
    ylabel("Amplitude")
    title("Signal généré sur la voie en phase");
figure;
    nexttile;   % Tracé signal x_Q
    plot(t, x_Q);
    xlim([0 t(length(t))]);
    xlabel("Temps (s)");
    ylabel("Amplitude")
    title("Signal généré sur la voie en quadrature");
figure;
    nexttile;   % Tracé signal x
    plot(t, x);
    xlim([0 t(length(t))]); 
    xlabel("Temps (s)");
    ylabel("Amplitude")
    title("Signal transmis sur fréquence porteuse");

figure;
    nexttile;   % Tracé DSP x_I
    semilogy (f, fftshift(DSP_x_I));
    xlabel("Fréquence (Hz)");
    ylabel("DSP");
    title("DSP du signal généré sur la voie de phase");
figure;
    nexttile;   % Tracé DSP x_Q
    semilogy (f, fftshift(DSP_x_Q));
    xlabel("Fréquence (Hz)");
    ylabel("DSP");
    title("DSP du signal généré sur la voie de quadrature");
figure;
    nexttile;   % Tracé DSP x
    semilogy(f, fftshift(DSP_x));
    xlabel("Fréquence (Hz)");
    ylabel("DSP");
    title("DSP du signal transmis sur fréquence porteuse");

figure;
    nexttile;   % Tracé TEB
    semilogy(EbsurN0, teb);
    xlabel("Eb/N0 (dB)");
    ylabel("TEB");
    title("TEB en fonction du rapport signal à bruit par bit");

figure;
    nexttile;   % Tracé TEB simulé vs théorique
    semilogy(EbsurN0, teb);
    hold on;
    semilogy(EbsurN0, teb_th);
    hold off;
    xlabel("Eb/N0 (dB)");
    ylabel("TEB");
    legend("TEB simulé", "TEB théorique");
    title("Comparaison des TEB simulé et théorique");
