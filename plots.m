clear all, close all, clc
Jmin = 10;
dt = 0.01;% time step
maxpop = 25;%  population size
maxgen = 12; %number of generations
s = tf('s');
G = 1/(s*(s*s+s+1));
rng(1,'twister');
load randpop.mat
options = optimoptions(@ga, 'PopulationSize', maxpop, 'MaxGenerations', maxgen, 'InitialPopulation', population, 'OutputFcn', @outfun);
[x, fval, exitflag, output, population, score] = ga(@(K)pidcalc(G, dt, K), 3, -eye(3), zeros(3,1), [], [], [], [], [], options); 
%%
load history.mat
for k=1:maxgen %sorted all columns which contain scores
   sortedcost(:,k) = sort(cost(:,k)); 
end
imagesc(log(sortedcost(:,1:maxgen)));
set(gcf, 'Position', [100 100 600 300]);
set(gcf, 'PaperPositionMode', 'auto');
print('-deps2', '-loose', '../../figures/GAPID1');
%%
figure 
hold on
    for k=1:maxgen
        for j=1:maxpop
              scatter3(history(j,1,k), history(j,2,k), history(j,3,k), 15, [(maxgen-j)/maxgen 0.25 maxgen/gen], 'filled');               
        end
    end
    [B, I] = sort(cost(:, 1:maxgen));
    scatter3(history(I(1), 1, maxgen), history(I(2), 2, maxgen), history(I(3), 3, maxgen), 100, [0 0 0], 'filled')
      view(69, 24)
      box on
      xlabel('P')
      ylabel('I')
      zlabel('D')
set(gcf,'Position',[100 100 350 250])
set(gcf,'PaperPositionMode','auto')
print('-deps2', '-loose', '../../figures/GAPID2');
%%
% we now compare step responses for the first and last generation 
%Generation 1
gen = 1;
t = 0:dt:20;
for k=1:maxpop
     K = history(k, 1, gen)+history(k, 2, gen)*s + history(k, 3, gen)*s/(1+0.001*s);%approximates derivative to avoid instability in step response causd by large increase due to derivative
     L = series(K,G); 
     cl = feedback(L, 1);%unity feedback
    [y, t] = step(cl, t);
    plot(t, y, 'LineWidth', 1.2);
end
box on, grid on
set(gcf,'Position',[100 100 550 250])
set(gcf,'PaperPositionMode','auto')
print('-deps2', '-loose', '../../figures/GAPID3');
%%
% Generation 10
gen = 10;
t = 0:dt:20;
for k=1:maxpop
     K = history(k, 1, gen)+history(k, 2, gen)*s + history(k, 3, gen)*s/(1+0.001*s);%approximates derivative to avoid instability in step response causd by large increase due to derivative
     L = series(K,G); 
     cl = feedback(L, 1);%unity feedback
    [y, t] = step(cl, t);
    plot(t, y, 'LineWidth', 1.2);
end
box on, grid on
set(gcf,'Position',[100 100 550 250])
set(gcf,'PaperPositionMode','auto')
print('-deps2', '-loose', '../../figures/GAPID4');
%%
for g=1:maxgen
     K = history(I(1), 1, g)+history(I(2), 2, g)*s + history(I(3), 3, g)*s/(1+0.001*s);%approximates derivative to avoid instability in step response causd by large increase due to derivative
     L = series(K,G);
     cl = feedback(L, 1);
     [y, t] = step(cl, t);
     subplot(3,1,1), hold on
     plot(t, y, 'LineWidth', 1+0.1*g, 'Color', [(maxgen-g)/maxgen 0 g/maxgen], 'LineWidth', 1.2);
     box on; grid on
     subplot(3,1,2), hold on
     cltf = K/(1+K*G);
     [u, t] = lsim(K, 1-y, t);
     plot(t, u, 'LineWidth', 1+0.1*g, 'Color', [(maxgen-g)/maxgen 0 g/maxgen], 'LineWidth', 1.2);
     q = 1;
     r = 0.001;
     J = dt*cumsum(q*(1-y(:)).^2 + r*u(:).^2);
     if (J<Jmin)
        Jmin = J; 
     end
end
disp(Jmin)
box on, grid on
set(gcf,'Position',[100 100 550 350])
set(gcf,'PaperPositionMode','auto')
print('-deps2', '-loose', '../../figures/GAPID5');