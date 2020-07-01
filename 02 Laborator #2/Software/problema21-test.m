%% 1

% 0.0

%% 2

clc;
noise;

% simuland filtrul de ordinul II putem observa ca spectrul acestuia se
% comporta ca un filtru trece-jos daca partea reala a polilor este mai mare
% decat partea reala a zeroului

%  atunci cand partea reala a polilor este mai mica decat partea reala a
%  zeroului, atunci spectrul filtrului incepe s? se comporte ca un filtru
%  trece-sus

% evident, in conditiile in care polii parasesc cercul unitar,
% sistemul devine instabil

% de asemenea se poate observa ca departarea polilor unul de celalalt
% cauzeaza un varf de rezonanta in spectrul filtrului

%% 3

% pentru a obtine un filtru trece-jos este necesar ca partea reala a 
%polilor sa fie cel putin mai mare(preferabil mult mai mare) decat partea 
%reala a zeroului. 

% de asemenea, pentru a obtine o comportare mai buna este de preferat ca 
% polii sa fie cat pai apropiati unul de celalalt(ideal ar fi ca partea 
% lor imaginara sa fie 0).

%% 4

% pentru a putea obtine un varf de rezonanta in punctul w = 1, partea reala
% a polilor trebuie sa fie cat mai aproape de 0.5, iar partea imaginara cat
% mai aproape de extremitatile cercului unitar

%% 5

% pozitia zeroului nu pare a afecta comportamentul sistemului. Mai 
% degraba este imortanta pozitia polilor fata de zerou.