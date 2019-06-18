%% Mjerenje signala - testni primjer

f = 5e6;                        % frekvencija 5 MHz
T = 1/f;                        % period
kut = 2*pi/3;                   % fazni pomak

t = 0:1e-8:1e-6 - 1e-8;
x = sin(2*pi*f*t + kut);
x = awgn(x, 20);                % SNR = 20dB
plot(t, x);
xlabel('Vremenski osdjecak 1us');

comb = zeros(1, 100);
for i=1:10:100
    comb(i) = x(i);
end
hold on;
stem(t, comb);
ylabel('x(t), x(n)');
title('Signal i deset uzoraka na odsjecku');
hold off;

%% Mjerenje signala - 15 signala razlicitih faznih pomaka za svaku od 15 frekvencija

f = linspace(3e6, 5e6, 15);            % 15 frekvencija od 3 MHz do 5 MHz
faza = linspace(0, 2*pi, 15);           % 15 faznih pomaka od 0 rad do 2pi rad
t = 0:1e-8:1e-6 - 1e-8;

X = zeros(15,100,15);                   % frekvencija x vrijeme x faza
C = zeros(15,100,15);

for i = 1:size(C, 1)                            % iterator frekvencije
    for j = 1:size(C, 3)                        % iterator faze
        X(i,:,j) = awgn(sin(2*pi*f(i)*t + faza(j)), 20);
        for k = 1:10:100
            C(i,k,j) = X(i,k,j);
        end
    end
end

%% Mjerenje signala - prikaz rezultata mjerenja 15x15 (225) signala

for i = 1:size(C, 1)
    for j = 1:size(C, 3)
        plot(t, X(i,:,j));
        hold on
        stem(t, C(i,:,j));
        hold off
        ylabel('x(t), x(n)');
        xlabel('Vremenski osdjecak 1us');
        title('Signal i deset uzoraka na odsjecku');
        pause;                                                   % iteracija se nastavlja s ENTER, prekida s CTRL+C
        msg = sprintf('Frekvencija: %f MHz | Faza: %f rad', f(i)/1e6, faza(j));
        disp(msg);
    end
end

%% Kovarijancijska matrica
% Za svaki skup signala jedne frekvencije izracunaj mean (mu) i kov. matricu (sigma)

for i = 1:size(C,1)                % iterator frekvencije

    sigma(:,:,i) = zeros(length(t), length(t));
    mu(:,i) = mean(C(i,:,:), 3);
    
    for k = 1:size(C, 3)              % iterator faze
        
       sigma(:,:,i) = sigma(:,:,i) + ((((C(i,:,k)-mu(:,i))))*(C(i,:,k)-mu(:,i))');
        
    end
    
    sigma(:,:,i) = sigma(:,:,i)./size(C, 3);
    
end

%% PCA, svojstveni vektori/vrijednosti

for i = 1:size(sigma, 3)
    
    [V(:,:,i),D(:,:,i)] = eig(sigma(:,:,i));
    
end

%% Zajednicka kovarijancijska matrica

sigma_un = (1./size(sigma, 3))*sum(sigma,3);
sigma_un = (1/2)*(sigma_un + sigma_un');
[V_un, D_un] = eig(sigma_un);

%% Mjerna matrica

M = 10;
phi = (V_un(:, end-M+1:end)');

%% Rekonstrukcija
plot_stem = 1;

for i = 1:size(C, 1)                    % iterator frekvencije
    for k = 1:size(C, 3)                % iterator faze
       
        C_test = C(i,:,k)';
        y = phi*C_test;          
        C_rec(:,:,k) = phi'*y;          % phi*phi' ~ I
        alfa_rec(:,:,k) = V(:,:,k)'*C_rec(:,:,k);
        sig_rec = V(:,:,k)* alfa_rec(:,:,k);
        
        if plot_stem
            stem(t, C_test, 'b');
            hold on;
            stem(t, real(sig_rec), 'r');
            hold off;
            xlabel('Vremenski osdjecak 1us');
            title('Signal i deset uzoraka na odsjecku');
            legend('Izvorni signal', 'Rekonstruirani signal')
            pause;                          % iteracija se nastavlja s ENTER, prekida s CTRL+C
        else
            plot(t, C_test, 'b');
            hold on;
            plot(t, real(sig_rec), 'r');
            hold off;
            xlabel('Vremenski osdjecak 1us');
            title('Signal i deset uzoraka na odsjecku');
            legend('Izvorni signal', 'Rekonstruirani signal')
            pause;                          % iteracija se nastavlja s ENTER, prekida s CTRL+C
        end
        
    end
end









