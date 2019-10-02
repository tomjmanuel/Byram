%%
x = -1000:1000;
f0 = 700;
y = sin(2*pi.*x./f0);
%y = padarray(y,[50 1]);

G = gausswin(max(size(x))-100);
G = G -min(G(:));
G = padarray(G,50);
close all
sig = G.*y';
sigp = padarray(sig,10000);
plot(sigp)

%%
figure
S = fft(sigp);
plot(log10(abs(S)))