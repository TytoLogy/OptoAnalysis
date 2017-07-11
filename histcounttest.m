S = 5:5:95;
edges = 0:10:100;

[N, E] = histcounts(S, edges);

figure(10)
centers = edges(1:end-1) + diff(edges)/2;
bar(centers, N, 1)


figure(11)
psth(S, 10, [0 100])