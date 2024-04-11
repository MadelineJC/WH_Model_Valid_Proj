%% Antia et al. (1994)

% Declaring the variables of the system
syms P I r k p o

% Setting up equations
dP = r*P-k*P*I
dI = p*I*(P/(P + o))

% Calculating Jacobian matrix
J = jacobian([dP, dI], [P, I])

% [                           r - I*k,          -P*k]
% [ (I*p)/(P + o) - (I*P*p)/(P + o)^2, (P*p)/(P + o)]

%% Fenton and Perkins (2010)

% Declaring the variables of the system
syms P I r B h e d

% Setting up equations
dP1 = r*P - I*B*P
dI1 = e*I*B*P - d*I

dP2 = r*P - I*(B*P/(1 + h*B*P))
dI2 = e*I*(B*P/(1 + h*B*P)) - d*I

dP3 = r*P - I*(B*P^2/(1 + h*B*P^2))
dI3 = e*I*(B*P^2/(1 + h*B*P^2)) - d*I

% Looking for equilibria

%% Mod. Fenton and Perkins (2010)
%% Type I

% Declaring the variables of the system
syms P I r B e d b

% Setting up equations
dP = r*P - I*B*P
dI = b + e*I*B*P - d*I

% Looking for equilibria
dP_0P = solve(dP, P) % P = 0
% dP_0I = solve(dP, I) % I = r/B
% dI_0P = solve(dI, P) % P = -(b - I*d)/(B*I*e)
dI_0I = solve(dI, I) % I = b/(d - B*P*e)

dI_New = subs(dI, P, 0)
dI_Eq = solve(dI_New, I) % I* = b/d

dP_New = subs(dP, I, dI_Eq)
dP_Eq = solve(dP_New, P) % P* = 0


% Calculating Jacobian matrix
J = jacobian([dP, dI], [P, I])

% [ r - B*I,      -B*P]
% [   B*I*e, B*P*e - d]

d1dP = r - B*I
d1dI = -B*P
d2dP = B*I*e
d2dI = B*P*e - d

JSub = subs([d1dP, d1dI, d2dP, d2dI],[P, I], [0, b/d])
JSub_Eig = eig(JSub)

%% Type II

% Declaring the variables of the system
syms P I r B h e d b

% Setting up equations
dP = r*P - I*(B*P/(1 + h*B*P))
dI = b + e*I*(B*P/(1 + h*B*P)) - d*I

% Looking for equilibria
dP_0P = solve(dP, P) % P = -(r - B*I)/(B*h*r)
dP_0P = -(r - B*I)/(B*h*r)
dI_0I = solve(dI, I) % I = b/(d - (B*P*e)/(B*P*h + 1))

dI_New = subs(dI, P, dP_0P)
dI_Eq = solve(dI_New, I) % I* = (e*r - B*b*h)/(B*e - B*d*h)

dP_New = subs(dP, I, dI_Eq)
dP_Eq = solve(dP_New, P) % P* = -(B*b - d*r)/(B*e*r - B*d*h*r)
dP_Eq = -(B*b - d*r)/(B*e*r - B*d*h*r)


% Calculating Jacobian matrix
J = jacobian([dP, dI], [P, I])

% [ r - (B*I)/(B*P*h + 1) + (B^2*I*P*h)/(B*P*h + 1)^2,      -(B*P)/(B*P*h + 1)]
% [ (B*I*e)/(B*P*h + 1) - (B^2*I*P*e*h)/(B*P*h + 1)^2, (B*P*e)/(B*P*h + 1) - d]

d1dP = r - (B*I)/(B*P*h + 1) + (B^2*I*P*h)/(B*P*h + 1)^2
d1dI = -(B*P)/(B*P*h + 1)
d2dP = (B*I*e)/(B*P*h + 1) - (B^2*I*P*e*h)/(B*P*h + 1)^2
d2dI = (B*P*e)/(B*P*h + 1) - d

JSub = subs([d1dP, d1dI, d2dP, d2dI],[P, I], [dP_Eq, dI_Eq])
JSub_Eig = eig(JSub)
