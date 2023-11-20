clear all;
close all;

Fe = 24000;         % Fréquence d'échantillonnage
Rb = 3000;          % Débit binaire de la transmission
N_bits = 50000;     % Nombre de bits à transmettre
f_p = 2000;         % Fréquence porteuse
EbsurN0 = 0:1:6;    % Rapports signal sur bruit à tester
alpha = 0.35;       % Roll off
L = 10;             % Span

Te = 1 / Fe;                                % Période d'échantillonnage
Tb = 1 / Rb;                                % Période par bits
Ns = fix(Tb / Te);                          % Facteur de suréchantillonnage
bits = randi([0, 1], 1, N_bits);            % Bits à transmettre


%% Génération du signal
symboles = (bits == 1) - (bits == 0);                     % Mapping BPSK
mat_kron = [1 zeros(1, Ns-1)];
h = rcosdesign(alpha, L, Ns);
Nt = L / 2;
retard = Nt * Ns;

Mod = kron(symboles, mat_kron);                            % Suréchantillonnage
Mod_protection_retard = cat(2, Mod, zeros(1, retard));     % Gestion du retard
x = filter(h, 1, Mod_protection_retard);                   % Filtrage
x = x(1, retard+1 : end);                                  % Gestion du retard
phi=40*pi/180;

[~, ~, teb_exp_sc_0] = transmission_erreur_phase_bruit(bits, Ns,h, N_bits, x, retard, 0);
[teb_theo_40, teb_exp_ce_40, teb_exp_sc_40] = transmission_erreur_phase_bruit(bits, Ns,h, N_bits, x, retard, 40/180*pi);
[teb_theo_100, teb_exp_ce_100, teb_exp_sc_100] = transmission_erreur_phase_bruit(bits, Ns,h, N_bits, x, retard, 100/180*pi);


figure
semilogy(teb_theo_40);
hold on;
semilogy(teb_exp_sc_40);
hold off;
title("TEB  avec erreur de phase (Erreur de phase = 40°)");
xlabel('Eb/N0 (dB)');
ylabel('Taux d''erreur binaire');
legend('théorique', 'experimentale');

figure
semilogy(teb_exp_sc_40);
hold on;
semilogy(teb_exp_sc_100);
hold off;
title("TEB  avec erreur de phase (Erreur de phase = 40° et 100°)");
xlabel('Eb/N0 (dB)');
ylabel('Taux d''erreur binaire');
legend('erreur de 40° ', 'erreur de 100°');


figure
semilogy(teb_exp_sc_0);
hold on;
semilogy(teb_exp_sc_40);
hold off;
title("TEB  avec erreur de phase (Erreur de phase = 40° et sans erreur de phase)");
xlabel('Eb/N0 (dB)');
ylabel('Taux d''erreur binaire');
legend('sans erreur ', 'erreur de 40°');


figure
semilogy(teb_exp_sc_0);
hold on;
semilogy(teb_exp_sc_0);
semilogy(teb_exp_ce_40);
semilogy(teb_exp_sc_100);
semilogy(teb_exp_ce_100);

hold off;
title("TEB  avec correction (Erreur de phase = 40° et 100°)");
xlabel('Eb/N0 (dB)');
ylabel('Taux d''erreur binaire');
legend('sans erreur ', 'erreur de 40°','correction erreur de 40°', 'erreur de 100°','correction erreur de 100°');

[teb_theo_40, teb_exp_codage_40, teb_exp_sc_40] = transmission_erreur_codage(bits, Ns,h, N_bits, symboles, retard, 40/180*pi);
[teb_theo_100, teb_exp_codage_100, teb_exp_sc_100] = transmission_erreur_codage(bits, Ns,h, N_bits, symboles, retard, 100/180*pi);


figure
semilogy(teb_exp_codage_100);
hold on;
semilogy(teb_exp_sc_0);
hold off;
title("");
xlabel('Eb/N0 (dB)');
ylabel('Taux d''erreur binaire');
legend('codage 100 ', 'sans codage');



%% Enveloppe convexe associée au signal

function [teb_theorique,teb_correction , teb_sans_correction] = transmission_erreur_phase_bruit(bits, Ns, h, N_bits, x_e,retard , phi)

    EbsurN0 = 0:1:6;    % Rapports signal sur bruit à tester

    teb_sans_correction = zeros(1, length(EbsurN0));
    teb_correction = zeros(1, length(EbsurN0));
    teb_theorique = zeros(1,length(EbsurN0));
    for i=1:length(EbsurN0)
        %% Canal de propagation à bruit complexe
        M = 2;                                                  % Ordre de la modulation (BPSK)
        Px_e = mean(abs(x_e).^2);                               % Puissance du signal transmis
        sigmacarre = Px_e * Ns / (2 * log2(M) * 10.^(EbsurN0(i)/10)); % Puissance du bruit pour le rapport Eb/N0 souhaité
        bruit = sqrt(sigmacarre) * (randn(1, length(x_e)) + 1i * randn(1, length(x_e))); % Bruit complexe
        y = x_e + bruit;                                        % Signal bruité
       
       %Ajout de l'erreur de phase 
    
       y = y *exp(1i *phi);
    
    
    
        %% Filtre de réception
        y_protection_retard = cat(2, y, zeros(1, retard));  % Gestion du retard
        z = filter(h, 1, y_protection_retard);              % Retour en bande de base
        z = z(1, retard+1 : end);                           % Gestion du retard
        
        %% Échantillonnage
        n0 = 1;                             % Instant d'échantillonnage (retard géré avant)
        z_echantillonne = z(n0 : Ns : end); % Échantillonnage à n0+m*Ns
        
    
        %% Estimation de l'erreur phi chapeau et application de la correction
         phi_chapeau = 1/2 * angle(sum(z_echantillonne.^2));
         z_correction_erreur =z_echantillonne* exp(-1i*phi_chapeau);
    
           
        

        
        %% Décision et démapping

        %avec correction d'erreur
        bits_demap_ce = 1*(z_correction_erreur >= 0);  
        %sans correction d'erreur
        bits_demap = 1*(z_echantillonne >= 0);     


        
        %% Taux d'erreur binaire
        teb_sans_correction(i)=sum(bits_demap ~= bits) / N_bits;
        teb_correction (i) = sum(bits_demap_ce ~= bits) / N_bits;
        teb_theorique(i) = qfunc(sqrt(2*10.^(EbsurN0(i)/10)*cos(phi)^2));

    end 
end 

function [teb_theorique,teb_correction , teb_sans_correction] = transmission_erreur_codage(bits, Ns, h, N_bits, x_e,retard , phi)

     xcode(1)=x_e(1);

       for k = 2 : length(x_e)
            xcode(k)= xcode(k-1)*x_e(k);
       end 


    mat_kron = [1 zeros(1, Ns-1)];
    Mod = kron(xcode, mat_kron);                            % Suréchantillonnage
    Mod_protection_retard = cat(2, Mod, zeros(1, retard));     % Gestion du retard
    x_e = filter(h, 1, Mod_protection_retard);                   % Filtrage
    x_e = x_e(1, retard+1 : end);                                  % Gestion du retard

    EbsurN0 = 0:1:6;    % Rapports signal sur bruit à tester

    teb_sans_correction = zeros(1, length(EbsurN0));
    teb_correction = zeros(1, length(EbsurN0));
    teb_theorique = zeros(1,length(EbsurN0));
    for i=1:length(EbsurN0)
        %% Canal de propagation à bruit complexe
        M = 2;                                                  % Ordre de la modulation (BPSK)
        Px_e = mean(abs(x_e).^2);                               % Puissance du signal transmis
        sigmacarre = Px_e * Ns / (2 * log2(M) * 10.^(EbsurN0(i)/10)); % Puissance du bruit pour le rapport Eb/N0 souhaité
        bruit = sqrt(sigmacarre) * (randn(1, length(x_e)) + 1i * randn(1, length(x_e))); % Bruit complexe
        y = x_e + bruit;                                        % Signal bruité
       
       %Ajout de l'erreur de phase 
    
       y = y *exp(1i *phi);
    
    
    
        %% Filtre de réception
        y_protection_retard = cat(2, y, zeros(1, retard));  % Gestion du retard
        z = filter(h, 1, y_protection_retard);              % Retour en bande de base
        z = z(1, retard+1 : end);                           % Gestion du retard
        
        %% Échantillonnage
        n0 = 1;                             % Instant d'échantillonnage (retard géré avant)
        z_echantillonne = z(n0 : Ns : end); % Échantillonnage à n0+m*Ns
        
    
        %% Estimation de l'erreur phi chapeau et application de la correction
         phi_chapeau = 1/2 * angle(sum(z_echantillonne.^2));
         z_correction_erreur =z_echantillonne* exp(-1i*phi_chapeau);
    
       z_decode(1)=z_correction_erreur(1);  
       for k = 2 : length(z_correction_erreur)
            z_decode(k)= z_decode(k-1)*z_correction_erreur(k);
       end 


        
        %% Décision et démapping

        %avec correction d'erreur
        bits_demap_ce = 1*(z_decode >= 0);  
        %sans correction d'erreur
        bits_demap = 1*(z_echantillonne >= 0);     


        
        %% Taux d'erreur binaire
        teb_sans_correction(i)=sum(bits_demap ~= bits) / N_bits;
        teb_correction (i) = sum(bits_demap_ce ~= bits) / N_bits;
        teb_theorique(i) = qfunc(sqrt(2*10.^(EbsurN0(i)/10)*cos(phi)^2));

    end 
end 



