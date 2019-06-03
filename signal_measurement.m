%% Mjerenje signala - testni primjer

f = 5e6;                        % frekvencija 5 MHz
T = 1/f;                        % period
kut = 2*pi/3;                   % fazni pomak

t = 0:1e-8:1e-6;
x = sin(2*pi*f*t + kut);
x = awgn(x, 20);                % SNR = 20dB
plot(t, x);
xlabel('Vremenski osdjecak 1us');

comb = zeros(1, 101);
for i=1:10:101
    comb(i) = x(i);
end
hold on;
stem(t, comb);
ylabel('x(t), x(n)');
title('Signal i deset uzoraka na odsjecku');
hold off;

%% Mjerenje signala - 15 signala razlicitih faznih pomaka za svaku od 15 frekvencija

f = linspace(5e6, 10e6, 15);            % 15 frekvencija od 5 MHz do 10 MHz
faza = linspace(0, 2*pi, 15);           % 15 faznih pomaka od 0 rad do 2pi rad
t = 0:1e-8:1e-6;

X = zeros(15,101,15);                   % frekvencija x vrijeme x faza
C = zeros(15,101,15);

for i = 1:15                            % iterator frekvencije
    for j = 1:15                        % iterator faze
        X(i,:,j) = awgn(sin(2*pi*f(i)*t + faza(j)), 20);
        for k = 1:10:101
            C(i,k,j) = X(i,k,j);
        end
    end
end

%% Mjerenje signala - prikaz rezultata mjerenja 15x15 (225) signala

for i = 1:15
    for j = 1:15
        plot(t, X(i,:,j));
        hold on
        stem(t, C(i,:,j));
        hold off
        ylabel('x(t), x(n)');
        xlabel('Vremenski osdjecak 1us');
        title('Signal i deset uzoraka na odsjecku');
        pause;                                                   % iteracija se nastavlja s ENTER, prekida s CTRL+C
        msg = sprintf('Frekvencija: %f | Faza: %f', f(i), faza(j));
        disp(msg);
    end
end



