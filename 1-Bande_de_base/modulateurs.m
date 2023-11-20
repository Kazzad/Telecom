clear all;
close all;

Fe = 24000;     % Fréquence d'échantillonage
Rb = 3000;      % Débit binaire de la transmission
N_bits = 30;  % Nombre de bits à transmettre

Te = 1 / Fe;                        % Période d'échantillonage
Tb = 1 / Rb;                        % Période par bits
Ns = fix (Tb / Te);                 % Facteur de suréchantillonage
bits = randi ([0, 1], 1, N_bits);   % Bits à transmettre

%% 2. ÉTUDE DE MODULATEURS BANDE DE BASE

%% === Modulateur 1
% Constantes
M1 = 2; % Nb de symboles
Rs1 = Rb / log2 (M1);
Ns1 = Fe / Rs1; % Facteur de suréchantillonage
% Mapping
symboles_1 = (bits==1)-(bits==0); % Signal mappé binaire à moyenne nulle
% Suréchantillonage
mat_kron = zeros (1, Ns1);
mat_kron (1) = 1;
Mod_1 = kron (symboles_1, mat_kron); % Signal suréchantilloné
% Filtrage
h1 = ones (1, Ns1);
Signal_Mod_1 = filter (h1, 1, Mod_1); % Signal généré
DSP1 = pwelch (Signal_Mod_1, [], [], [], Fe, 'twosided'); % DSP estimée
f1 = linspace (-Fe/2, Fe/2, length (DSP1));
DSP1_theo = Tb * sinc (f1*Tb) .^ 2; % DSP théorique

figure ("Name", "Modulateur 1");
    % Signal généré
    nexttile;
    plot (Signal_Mod_1);
    xlim([0, length(Signal_Mod_1)]);
    ylim ([-1.2, 1.2]);
    xlabel ("Temps (s)");
    ylabel ("Signal");
    title ("Signal généré");
    % DSP estimée
    nexttile;
    semilogy (f1, fftshift(DSP1));
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    title ("DSP estimée");
    % DSP estimée et DSP théorique
    nexttile;
    semilogy (f1, fftshift(DSP1));
    hold on;
    semilogy (f1, DSP1_theo);
    hold off;
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    legend ("estimée", "théorique");
    title ("DSPs estimée et théorique");

%% === Modulateur 2
% Constantes
M2 = 4; % Nb de symboles
Rs2 = Rb / log2 (M2);
Tb2 = 1 / Rs2;
Ns2 = Fe / Rs2; % Facteur de suréchantillonage
% Mapping
bits_reshape = reshape (bits, 2, N_bits / 2)'; % Bits par paire sur 2 colonnes
symboles_2 = -3 * (bits_reshape(:, 1) == 0 & bits_reshape(:, 2) == 1);              %   01 -> -3
symboles_2 = symboles_2 + (bits_reshape(:, 1) == 0 & bits_reshape(:, 2) == 0);      %   00 ->  1
symboles_2 = symboles_2 + 3 * (bits_reshape(:, 1) == 1 & bits_reshape(:, 2) == 0);  %   10 ->  3
symboles_2 = symboles_2 - (bits_reshape(:, 1) == 1 & bits_reshape(:, 2) == 1);      %   11 -> -1
symboles_2 = symboles_2'; % Signal mappé 4-aires à moyenne nulle
% Suréchantillonage
mat_kron = zeros (1, Ns2);
mat_kron (1) = 1;
Mod_2 = kron (symboles_2, mat_kron);  % Signal suréchantilloné
% Filtrage
h2 = ones (1, Ns2);
Signal_Mod_2 = filter (h2, 1, Mod_2); % Signal généré
DSP2 = pwelch (Signal_Mod_2, [], [], [], Fe, 'twosided'); % DSP estimée
f2 = linspace (-Fe/2, Fe/2, length (DSP2));
DSP2_theo = Tb2 * sinc (f2*Tb2) .^ 2; % DSP théorique

figure ("Name", "Modulateur 2");
    % Signal généré
    nexttile;
    plot (Signal_Mod_2);
    xlim([0, length(Signal_Mod_2)]);
    ylim ([-3.5, 3.5]);
    xlabel ("Temps (s)");
    ylabel ("Signal");
    title ("Signal généré");
    % DSP estimée
    nexttile;
    semilogy (f2, fftshift(DSP2));
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    title ("DSP estimée");
    % DSP estimée et DSP théorique
    nexttile;
    semilogy (f2, fftshift(DSP2));
    hold on;
    semilogy (f2, DSP2_theo);
    hold off;
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    legend ("estimée", "théorique");
    title ("DSPs estimée et théorique");

%% === Modulateur 3
% (on réutilise les calculs déjà fait pour le modulateur 1 qui utilise le même mapping)
% Constantes
Ns3 = Ns1;
% Suréchantillonage
Mod_3 = Mod_1; % Signal suréchantilloné (le même que pour le mod 1)
alpha = 0.8; % Roll off
% Filtrage
h3 = rcosdesign(alpha, 10, Ns3);
Signal_Mod_3 = filter (h3, 1, Mod_3); % Signal généré
DSP3 = pwelch (Signal_Mod_3, [], [], [], Fe, 'twosided'); % DSP estimée
f3 = linspace (-Fe/2, Fe/2, length (DSP3));
sigma_carre = var (symboles_1);
DSP3_theo = zeros (1, length (DSP3));
DSP3_theo = DSP3_theo + (abs(f3) <= ((1-alpha)/(2*Tb))) * Tb;
DSP3_theo = DSP3_theo + Tb/2 *(1+cos (pi *Tb/alpha *(abs(f3)-(1-alpha)/(2*Tb)))) .* ((1-alpha)/(2*Tb) <= abs(f3) & abs(f3) <= (1+alpha)/(2*Tb));
DSP3_theo = DSP3_theo .* sigma_carre / Tb; % DSP théorique

figure ("Name", "Modulateur 3");
    % Signal généré
    nexttile;
    plot (Signal_Mod_3);
    xlim([0, length(Signal_Mod_3)]);
    ylim ([-0.7, 0.7]);
    xlabel ("Temps (s)");
    ylabel ("Signal");
    title ("Signal généré");
    % DSP estimée
    nexttile;
    semilogy (f3, fftshift(DSP3));
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    title ("DSP estimée");
    % DSP estimée et DSP théorique
    nexttile;
    semilogy (f3, fftshift(DSP3));
    hold on;
    semilogy (f3, DSP3_theo);
    hold off;
    xlabel ("Fréquence (Hz)");
    ylabel ("DSP");
    legend ("estimée", "théorique");
    title ("DSPs estimée et théorique");

figure ("Name", "DSP des 3 modulateurs");
    semilogy (f1, fftshift(DSP1));
    hold on;
    semilogy (f2, fftshift(DSP2));
    hold on;
    semilogy (f3, fftshift(DSP3));
    hold off;
    legend ("1", "2", "3");
    title ("DSP estimées des 3 modulateurs");

%% 3. ÉTUDE DES INTERFÉRENCES ENTRE SYMBOLE ET CRITÈRE DE NYQUIST

%% === Étude sans canal de propagation
hr1 = ones (1, Ns1);
sortie_filtre_sc_sb = filter (hr1, 1, Signal_Mod_1);

figure ("Name", "Sans canal");
    % Sortie du filtre
    nexttile;
    plot(sortie_filtre_sc_sb);
    ylim ([-8.5, 8.5]);
    xlabel ("Temps (s)");
    ylabel ("Signal");
    title ("Signal en sortie du filtre");
    % Réponse impulsionnelle
    nexttile;
    triangle = conv(h1,hr1);
    plot (triangle);
    xlabel ("Temps (s)");
    ylabel ("Réponse imp");
    title ("Réponse impulsionnelle de la chaîne de transmission");
    % Diagramme de l'oeil
    nexttile
    plot(reshape(sortie_filtre_sc_sb,Ns,length(sortie_filtre_sc_sb)/Ns));
    title ("Diagramme de l'oeil en sortie du filtre");

% 1. D'après réponse impulsionnelle et diagramme de l'oeil
n0 = Ns;
signal_demap = ones(1, N_bits);
m=8;
% Échantillonnage
for k = 0:(N_bits-1)
    signal_demap(k+1) = sortie_filtre_sc_sb(k*m+n0);
end
signal_demap = (signal_demap==8)+0*(signal_demap==-8);
taux_erreur_sc_sb = 1-sum(signal_demap ==bits)/length(signal_demap) % = 0

% 2. Modification de l'instant d'échantillonage 
critere_nyquist = 3;
signal_demap = ones(1, N_bits);
m=8;
for k = 0:(N_bits-1)
    signal_demap(k+1) = sortie_filtre_sc_sb(k*m+n0);
end 
signal_demap_biaise = (signal_demap==8)+0*(signal_demap==-8);
taux_erreur_sc_sb_biaise = 1-sum(signal_demap ==bits)/length(signal_demap) % = 1

%% === Étude avec canal de propagation sans bruit
% Canal de propagation à bande limitée BW = 8000Hz et BW = 1000Hz
[~, ~, ~, taux_erreur_ac_sb_8000] = canal_propagation_sans_bruit (8000,N_bits,bits,Ns,Fe,Signal_Mod_1,hr1,triangle); 
[~, ~, ~, taux_erreur_ac_sb_1000] = canal_propagation_sans_bruit (1000,N_bits,bits,Ns,Fe,Signal_Mod_1,hr1,triangle);
taux_erreur_ac_sb_8000 % = 0
taux_erreur_ac_sb_1000 % = 0.5 environ

function afficher_canal_propagation_sans_bruit (BW,Ns,Fe,triangle,rep_impulsionnelle,sortie_filtre_ac_sb,Filtre_PB)
    figure ("Name", "Avec canal, sans bruit (BW=" + BW + ")");
    % Réponse impulsionnelle
    nexttile;
    plot(rep_impulsionnelle);
    xlabel ("n");
    title ("Réponse impulsionnelle de la chaîne de transmission");
    % Diagramme de l'oeil
    nexttile;
    plot(reshape(sortie_filtre_ac_sb,Ns,length(sortie_filtre_ac_sb)/Ns));
    title ("Diagramme de l'oeil en sortie du filtre");
    % |H(f)H_r(f)| et |H_c(f)|
    nexttile;
    plot(linspace(-Fe,Fe,200),abs(fft(Filtre_PB,200)));
    hold on;
    plot(linspace(-Fe,Fe,200),abs(fftshift(fft(triangle,200))));
    xlabel ("Fréquence (Hz)");
    ylabel ("Amplitude");
    legend ("|H_c(f)|","|H(f)H_r(f)|");
    title("Réponses du filtre reception et mise en forme et du filtre canal");
end

function [rep_impulsionnelle, sortie_filtre_ac_sb, Filtre_PB, taux_erreur_ac_sb] = canal_propagation_sans_bruit (BW,N_bits,bits,Ns,Fe,Signal_Mod_1,hr1,triangle)
    % Filtre passe-bas
    f_coupure = 2*BW;
    f_coupure_tilde = f_coupure/Fe;
    ordre_filtre = 15;
    N = (ordre_filtre-1)/2;
    K = (-N:N)';
    Filtre_PB = 2*f_coupure_tilde*sinc(2*f_coupure_tilde*K);
    
    % Passage par le filtre passe-bas
    filtre_protection_retard = cat(2,Signal_Mod_1,zeros(1,N));
    filtre_protection_retard = filter(Filtre_PB',1,filtre_protection_retard);
    sortie_filtre_pb = filtre_protection_retard(1,N+1:end);
    
    % Passage par le filtre de réception
    sortie_filtre_ac_sb = filter (hr1, 1, sortie_filtre_pb);
    rep_impulsionnelle = conv(triangle,Filtre_PB);

    % Affichage
    afficher_canal_propagation_sans_bruit (BW,Ns,Fe,triangle,rep_impulsionnelle,sortie_filtre_ac_sb,Filtre_PB);

    % D'après réponse impulsionnelle et diagramme de l'oeil
    n0 = 8;
    signal_demap_ac_sb = ones(1, N_bits);
    epsilon = 0.8;
    m=8;
    for k = 0:(N_bits-1)
        signal_demap_ac_sb(k+1) = sortie_filtre_ac_sb(k*m+n0);
    end
    signal_demap_ac_sb = (signal_demap_ac_sb <=epsilon+8).*(signal_demap_ac_sb >=-epsilon+8)+0*(signal_demap_ac_sb==-8+epsilon).*(signal_demap_ac_sb==-8-epsilon);
    taux_erreur_ac_sb = 1-sum(signal_demap_ac_sb ==bits)/length(signal_demap_ac_sb);
end