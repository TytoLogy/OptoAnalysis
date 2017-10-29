
%  a = randi([50 100], [1000 1]); 
%  b = randi([-100 -50], [1000 1]);

muA = 10;
sdA = 5;
muB = 0;
sdB = sdA;
nI = 10000;
sampleN = 50;

a = normrnd(muA, sdA, nI, 1); 
b = normrnd(muB, sdB, nI, 1);

c = change_index_stat(a, b)';

figure(1)
subplot(311)
hist(a, 20);
subplot(312);
hist(b, 20);
subplot(313);
hist(c, 20);

% now boot strap
aa = normrnd(muA, sdA, sampleN, 1); 
bb = normrnd(muB, sdB, sampleN, 1);

D = bootChangeIndex(aa, bb, nI);

figure(2)
subplot(311)
hist(D.sampleA, 20);
subplot(312);
hist(D.sampleB, 20);
subplot(313);
hist(D.sampleC, 20);



