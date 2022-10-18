function J = pidcalc(G,dt,parms)

s = tf('s');
K = parms(1) + parms(2)/s + parms(3)*s/(1+.001*s);
Loop = series(K,G);
cl = feedback(Loop,1);
t = 0:dt:20;
[y,t] = step(cl,t);

cltf = K/(1+K*G);
u = lsim(K,1-y,t);

Q = 1;
R = .001;
J = dt*sum(Q*(1-y(:)).^2 + R*u(:).^2)

step(cl,t)
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
drawnow