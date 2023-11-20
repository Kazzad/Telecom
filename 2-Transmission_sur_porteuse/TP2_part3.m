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
    bits_pairs (i) = bits (2*i);
    bits_impairs (i) = bits (2*i-1);
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

%% Enveloppe convexe associée au signal
x_e = x_I + 1i*x_Q;

teb = zeros(1, length(EbsurN0));
for i=1:length(EbsurN0)
    %% Canal de propagation à bruit complexe
    M = 4;                                                  % Ordre de la modulation
    Px_e = mean(abs(x_e).^2);                               % Puissance du signal transmis
    sigmacarre = Px_e*Ns / (2*log2(M)*10.^(EbsurN0(i)/10)); % Puissance du bruit pour le rapport Eb/N0 souhaité
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
    z_I = real(z_echantillonne);        % Partie réelle
    z_Q = imag(z_echantillonne);        % Partie imaginaire

    %% Plot constellation
    scatterplot(z_echantillonne);
    xlabel("En phase");
    ylabel("En quadrature");
    title("Constellation en sortie de l'échantillonneur (Eb/N0=" + EbsurN0(i) + "dB)");
    
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

DSP_x_e = pwelch(x_e, [], [], [], Fe, 'twosided');                      % DSP estimée de x_e   
teb_th = 2*qfunc(sqrt(2*log2(M)*10.^(EbsurN0/10))*sin(pi/M))/log2(M);   % TEB théorique pour une chaine M-PSK
f = linspace(-Fe/2, Fe/2, length (DSP_x_e));                            % Échelle fréquentielle
constellation = symboles_impairs + 1i*symboles_pairs;

%% Figures
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
    nexttile;   % Tracé DSP x_e
    semilogy(f, fftshift(DSP_x_e));
    xlabel("Fréquence (Hz)");
    ylabel("DSP");
    title("DSP de l'enveloppe convexe associée au signal transmis sur fréquence porteuse");

% Constellation sortie de mapping
    scatterplot(constellation);
    xlim([-2 2]);
    ylim([-2 2]);
    xlabel("En phase");
    ylabel("En quadrature");
    title("Constellation en sortie du mapping");

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
    legend("TEB chaine p-b éq.", "TEB chaine transp. de f.");
    title("Comparaison des TEB des deux chaines");
