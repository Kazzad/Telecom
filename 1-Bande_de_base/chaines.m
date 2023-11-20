clear all;
close all;

bits = randi ([0, 1], 1, 10000);    % Bits à transmettre
EbsurN0_liste = linspace(0,8,9);    % Rapports signal à bruit par bit (de 1 à 8)
Fe = 24000;                         % Fréquence d'échantillonage
Rb = 3000;                          % Débit binaire de la transmission

afficher_plots_chaine_individuelles = 1;
afficher_plots_comparaisons_chaines = 1;

[Ns1, sortie_fr_sans_bruit1, taux_erreur_chaine_sb1, sortie_fr_avec_bruit1, taux_erreur_chaine_ab_est1, taux_erreur_chaine_ab_th1] = chaine (1, bits, EbsurN0_liste, Fe, Rb, afficher_plots_chaine_individuelles);
[Ns2, sortie_fr_sans_bruit2, taux_erreur_chaine_sb2, sortie_fr_avec_bruit2, taux_erreur_chaine_ab_est2, taux_erreur_chaine_ab_th2] = chaine (2, bits, EbsurN0_liste, Fe, Rb, afficher_plots_chaine_individuelles);
[Ns3, sortie_fr_sans_bruit3, taux_erreur_chaine_sb3, sortie_fr_avec_bruit3, taux_erreur_chaine_ab_est3, taux_erreur_chaine_ab_th3] = chaine (3, bits, EbsurN0_liste, Fe, Rb, afficher_plots_chaine_individuelles);

% Affichages de comparaison entre les chaines
if (afficher_plots_comparaisons_chaines == 1)
    % Ch1 vs ch2
    figure ("Name", "Comparaison TEB 1 et 2");
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_est1);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_th1);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_est2);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_th2);
        hold off;
        title ("TEB des chaines 1 et 2");
        xlabel ("E_b/N_0");
        ylabel ("TEB");
        legend ("estimé ch1", "théorique ch1", "estimé ch2", "théorique ch2");
    % Ch1 vs ch3
        figure ("Name", "Comparaison TEB 1 et 3");
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_est1);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_th1);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_est3);
        hold on;
        semilogy (EbsurN0_liste, taux_erreur_chaine_ab_th3);
        hold off;
        title ("TEB des chaines 1 et 3");
        xlabel ("E_b/N_0");
        ylabel ("TEB");
        legend ("estimé ch1", "théorique ch1", "estimé ch3", "théorique ch3");
end

% En sortie :
% Ns = facteur de suréchantillonnage
% (sans bruit)
% sortie_fr_sans_bruit = Sortie FdR
% taux_erreur_chaine_sb = TEB estimé
% (avec bruit)
% sortie_fr_avec_bruit = Sortie FdR. Chaque ligne correspond à une valeur de Eb/N0
% taux_erreur_chaine_ab_est/th = TEB estimé et théorique. Chaque valeur correspond à une valeur de Eb/N0
function [Ns, sortie_fr_sans_bruit, taux_erreur_chaine_sb, sortie_fr_avec_bruit, taux_erreur_chaine_ab_est, taux_erreur_chaine_ab_th] = chaine (n_chaine, bits, EbsurN0_liste, Fe, Rb, afficher_plots_chaine_individuelles)
    N_bits = length(bits);  % Nombre de bits à transmettre

    %% Modulateur
    M = nb_symboles (n_chaine); % Nb de symboles
    Rs = Rb / log2 (M);
    Ns = Fe / Rs; % Facteur de suréchantillonage
    Signal_Mod = signal_mis_en_forme (n_chaine, Ns, bits);
    
    %% Étude sans bruit
    % === Sortie du filtre de réception et signal démappé
    [sortie_fr_sans_bruit, demap_sans_bruit] = reception_demapping (n_chaine, Signal_Mod, Ns, N_bits);

    % === Taux d'erreur
    taux_erreur_chaine_sb = 1 - sum(demap_sans_bruit == bits) / length(demap_sans_bruit);
    
    %% Étude avec bruit 
    Px = mean(abs(Signal_Mod).^2); % Puissance du signal à bruiter
    
    taux_erreur_chaine_ab_est = zeros(1, length(EbsurN0_liste));
    taux_erreur_chaine_ab_th = zeros(1, length(EbsurN0_liste));
    sortie_fr_avec_bruit = zeros (length(EbsurN0_liste), length(Signal_Mod));
    
    for i = 1:length(EbsurN0_liste)
        EbsurN0 = EbsurN0_liste (i);

        % === Canal de propagation à bruit additif et gaussien
        sigmacarre = Px*Ns / (2*log2(M)*10^(EbsurN0/10));           % Calcul de la puissance du bruit sigma²
        bruit = sqrt(sigmacarre) * randn(1, length(Signal_Mod));    % Bruit ajouté 
        Signal_bruite = Signal_Mod + bruit;                         % Signal bruité

        % === Sortie du filtre de réception et signal démappé
        [sortie_fr_avec_bruit(i, :), demap_avec_bruit] = reception_demapping (n_chaine, Signal_bruite, Ns, N_bits);

        % === Taux d'erreur
        taux_erreur_chaine_ab_est(i) = 1 - sum(demap_avec_bruit == bits) / length(demap_avec_bruit);
        taux_erreur_chaine_ab_th(i) = teb_th (n_chaine, M, EbsurN0);
    end

    %% Affichage
    if (afficher_plots_chaine_individuelles == 1)
        taux_erreur_chaine_sb

        figure ("Name", "Chaine " + n_chaine + ": diagramme oeil (sb)");
                nexttile;
                plot (reshape (sortie_fr_sans_bruit, Ns, length (sortie_fr_sans_bruit) / Ns));
                xlabel ("Temps normalisé (t/N_s)");
                ylabel ("Amplitude");
                title ("Diagramme de l'oeil en sortie du FdR (ch" + n_chaine + ", sans bruit)");

        figure ("Name", "Chaine " + n_chaine + ": diagramme oeil (ab)");
            for i = 1:length(EbsurN0_liste)
                EbsurN0 = EbsurN0_liste (i);
                nexttile;
                plot (reshape (sortie_fr_avec_bruit(i, :), Ns, length(sortie_fr_avec_bruit(i, :)) / Ns));
                xlabel ("Temps normalisé (t/N_s)");
                ylabel ("Amplitude");
                title ("Diagramme de l'oeil en sortie du FdR (ch" + n_chaine + ", E_b/N_0=" + EbsurN0 + ")");
            end
    
        figure ("Name", "Chaine " + n_chaine + ": taux d'erreur (ab)");
            semilogy (EbsurN0_liste, taux_erreur_chaine_ab_est);
            hold on;
            semilogy (EbsurN0_liste, taux_erreur_chaine_ab_th);
            hold off
            title ("Taux d'erreur en fonction de E_b/N_0 (ch" + n_chaine + ")");
            legend ("estimé", "théorique");
    end
end


% Signal en sortie de filtre de réception et signal demappé
function [sortie_filtre_reception, signal_demap] = reception_demapping (n_chaine, entree_filtre_reception, Ns, N_bits)  
    % === Filtre de réception
    hr = filtre_reception (n_chaine, Ns);

    % === Filtrage
    sortie_filtre_reception = filter (hr, 1, entree_filtre_reception);
    
    % === Échantillonnage aux instants optimaux
    n0 = instants_optimaux (n_chaine, Ns);

    % === Démapping
    signal_demap = demapping (n_chaine, sortie_filtre_reception, n0, N_bits);
end

% Signal en sortie du filtre de mise en forme
function [Signal_Mod] = signal_mis_en_forme (n_chaine, Ns, bits)
    % === Mapping
    symboles = mapping (n_chaine, bits);
    
    % === Suréchantillonage
    mat_kron = zeros (1, Ns);
    mat_kron (1) = 1;
    Mod = kron (symboles, mat_kron); % Signal suréchantilloné

    % === Filtre de mise en forme
    h = ones (1, Ns); % Filtre

    % === Filtrage
    Signal_Mod = filter (h, 1, Mod);
end



% Filtre de réception
function [hr] = filtre_reception (n_chaine, Ns)
    if (n_chaine == 1) || (n_chaine == 3)
        hr = ones (1, Ns);
    else
        hr = ones (1, Ns / 2);
    end
end

% Mapping
function [symboles] = mapping (n_chaine, bits)
    if (n_chaine == 1) || (n_chaine == 2)
        % Mapping binaire à moyenne nulle
        symboles = (bits==1)-(bits==0);
    else
        % Mapping 4-aire à moyenne nulle
        bits_reshape = reshape (bits, 2, length(bits) / 2)';
        symboles = 3*(bits_reshape(:, 1) == 1 & bits_reshape(:, 2) == 0);               % 10 -> 3
        symboles = symboles + (bits_reshape(:, 1) == 0 & bits_reshape(:, 2) == 0);      % 00 -> 1
        symboles = symboles - (bits_reshape(:, 1) == 0 & bits_reshape(:, 2) == 1);      % 01 -> -1
        symboles = symboles - 3*(bits_reshape(:, 1) == 1 & bits_reshape(:, 2) == 1);    % 11 -> -3
        symboles = symboles';
    end
end

% Signal demappé
function [signal_demap] = demapping (n_chaine, signal, n0, N_bits)
    signal_demap = ones (1, N_bits); % taille de signalMod1 / n0
    if (n_chaine == 1) || (n_chaine == 2)
           m=8;
        % Échantillonnage
        for k = 0:(N_bits-1)
            signal_demap(k+1) = signal(k*m+n0);
        end
        % Décision
        signal_demap = 1*(signal_demap >= 0) + 0*(signal_demap < 0);
    else
            m=16;
        for k = 0:N_bits/2-1
            % Échantillonnage
            symb_k = signal(k*m+n0);
            % Décision
            % 1er bit
            signal_demap(2*k+1) = 1*(symb_k > 32);                                            % 1 si dans ]+inf,32[
            signal_demap(2*k+1) = signal_demap(2*k+1) + 0*((32 >= symb_k) && (symb_k > 0));   % 0 si dans [32,0[
            signal_demap(2*k+1) = signal_demap(2*k+1) + 0*((0 >= symb_k) && (symb_k > -32));  % 0 si dans [0,-32[
            signal_demap(2*k+1) = signal_demap(2*k+1) + 1*(-32 >= symb_k);                    % 1 si dans [-32,-inf[
            % 2nd bit
            signal_demap(2*k+2) = 0*(symb_k > 32);                                              % 0 si dans ]+inf,32[
            signal_demap(2*k+2) = signal_demap(2*k+2) + 0*((32 >= symb_k) && (symb_k > 0));       % 0 si dans [32,0[
            signal_demap(2*k+2) = signal_demap(2*k+2) + 1*((0 >= symb_k) && (symb_k > -32));      % 1 si dans [0,-32[
            signal_demap(2*k+2) = signal_demap(2*k+2) + 1*(-32 >= symb_k);                        % 1 si dans [-32,-inf[
        end
    end
end

function [taux_erreur_chaine_ab_th] = teb_th (n_chaine, M, EbsurN0)
    if (n_chaine == 2)
        taux_erreur_chaine_ab_th = qfunc(sqrt(10^(EbsurN0/10)));
    else
        taux_erreur_chaine_ab_th = 2*(M-1) / (M*log2(M)) * qfunc(sqrt(1/(M^2-1) * 6*log2(M) * (10^(EbsurN0/10))));
    end
end

% Nb de symboles utilisés
function [M] = nb_symboles (n_chaine)
    if (n_chaine == 1) || (n_chaine == 2)
        M = 2; % binaire
    else
        M = 4; % 4-aire
    end
end

function [n0] = instants_optimaux (n_chaine, Ns)
    if (n_chaine == 1)
        n0= Ns;
    else 
        if (n_chaine == 2)
            n0 = Ns/2;
        else
            n0 =Ns;
        end
    end
end
