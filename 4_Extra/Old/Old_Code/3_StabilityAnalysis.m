%% EEB1455 Project Script
%% THIS IS A MATLAB SCRIPT

% Declaring the variables of the system
syms P H r O h b c u

%% Setting up equations

E1_Prey = P.*r - H.*(O)
E1_Pred = b + H.*(c.*(O.*P) - u)

E2_Prey = P.*r - H.*(O.*P/(1 + h.*O.*P))
E2_Pred = b + H.*(c.*(O.*P/(1 + h.*O.*P)) - u)

E3_Prey = P.*r - H.*(O.*P^2/(1 + h.*O.*P^2))
E3_Pred = b + H.*(c.*(O.*P^2/(1 + h.*O.*P^2)) - u)

%% Calculating the Jacobian matrix
% Taking partial derivatives of each set of equations, which respect to
% each of the state variables 

J1 = jacobian([E1_Pred,E1_Prey],[P,H])
J2 = jacobian([E2_Pred,E2_Prey],[P,H])
J3 = jacobian([E3_Pred,E3_Prey],[P,H])

%% Equilibrium states and Stability Analyses %%

%% Type 1

% First, invasion equilibrium where parasite abundance = 0.

P1Equil_Invas = 0 % P = 0
E1_Pred_New = b + H.*(c.*(O.*0) - u)
H1Equil_Invas = solve(E1_Pred_New, H) % H = b/u

% Second, coexistence/infection equilibrium where parasite and host both
% exist in the system 

E1_Prey_New2 = solve(E1_Prey,P) % dP/dt = 0 when (1) P = 0 or (2) r = HO (H = r/O)
H1Equil_Coex = r/O
E1_Pred_New2 = subs(E1_Pred,H,HEquil_Coex)
P1Equil_Coex = solve(E1_Pred_New2,P) % P = -(b - (r*u)/O)/(c*r)

% Jacobian matricies & dominant eigenvalues
d1PreydP = H.*O.*c
d1PreydH = O.*P.*c - u
d1PreddP = r
d1PreddH = -O

% Invasion
J1Inv = subs([d1PreydP,d1PreydH,d1PreddP,d1PreddH],[P,H],[0,b/u])
% [ (O*b*c)/u, -u, r, -O]
J1InvEig = eig(J1Inv) % = (O*b*c)/u
% Unstable when (O*b*c)/u > 0

% Coexistence
J1Coex = subs([d1PreydP,d1PreydH,d1PreddP,d1PreddH],[P,H],[-(b - (r*u)/O)/(c*r),r/O])
% [ c*r, - u - (O*(b - (r*u)/O))/r, r, -O], but you can simplify the second
% term to (-bO/r)
J1CoexEig = r
% Stable when r < 0; BUT DOES THIS MAKE SENSE? MatLab says domeig is cr...

%% Type II

% First, invasion equilibrium where parasite abundance = 0.

P2Equil_Invas = 0 % P = 0
E2_Pred_New = b + H.*(c.*(O.*0) - u)
H2Equil_Invas = solve(E2_Pred_New, H) % H = b/u

% Second, coexistence/infection equilibrium where parasite and host both
% exist in the system 

E2_Prey_New2 = solve(E2_Prey,P) % P = -(r - H*O)/(O*h*r)
E2_Pred_New2 = subs([E2_Pred],[P],[-(r - H*O)/(O*h*r)])
H2Equil_Coex = solve(E2_Pred_New2,H) % H = (c*r - O*b*h)/(O*c - O*h*u)
P2_Prey_New2 = subs(E2_Prey,H,(c*r - O*b*h)/(O*c - O*h*u))
P2Equil_Coex = solve(P2_Prey_New2,P) % P = -(O*b - r*u)/(O*c*r - O*h*r*u)

% Jacobian matricies & dominant eigenvalues 
d2PreydP = H*((O*c)/(O*P*h + 1) - (O^2*P*c*h)/(O*P*h + 1)^2)
d2PreydH = (O*P*c)/(O*P*h + 1) - u
d2PreddP = r - (H*O)/(O*P*h + 1) + (H*O^2*P*h)/(O*P*h + 1)^2
d2PreddH = -(O*P)/(O*P*h + 1)

% Invasion
J2Inv = subs([d2PreydP,d2PreydH,d2PreddP,d2PreddH],[P,H],[0,b/u])
% [ (O*b*c)/u, -u, r - (O*b)/u, 0]
J2InvEig = eig(J2Inv) % = (O*b*c)/u
% Unstable when (O*b*c)/u > 0

% Coexistence
J2Coex = subs([d2PreydP,d2PreydH,d2PreddP,d2PreddH],[P,H],[-(O*b - r*u)/(O*c*r - O*h*r*u),(c*r - O*b*h)/(O*c - O*h*u)])
J2CoexEig = eig(J2Coex) % This is a mess, but when is is not; might have to just do this numerically, using Type I as a guide

%% Type III

% First, invasion equilibrium where parasite abundance = 0.

P3Equil_Invas = 0 % P = 0
E3_Pred_New = subs(E3_Pred,P,0)
H3Equil_Invas = solve(E3_Pred_New,H) % H = b/u

% Second, coexistence/infection equilibrium where parasite and host both
% exist in the system 
E3_Prey_New2 = solve(E3_Prey,P) % P = ((O*(O*H^2 - 4*h*r^2))^(1/2) + H*O)/(2*O*h*r) **OR** P = -((O*(O*H^2 - 4*h*r^2))^(1/2) - H*O)/(2*O*h*r)
E3_Pred_New2_A = subs([E3_Pred],[P],[((O*(O*H^2 - 4*h*r^2))^(1/2) + H*O)/(2*O*h*r)])
E3_Pred_New2_B = subs([E3_Pred],[P],[-((O*(O*H^2 - 4*h*r^2))^(1/2) - H*O)/(2*O*h*r)])
H3Equil_Coex_A = solve(E3_Pred_New2_A,H) % (c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) + O*b*c - 2*O*b*h*u)/(2*O*u*(c - h*u)) **OR** -(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) - O*b*c + 2*O*b*h*u)/(2*O*u*(c - h*u))
H3Equil_Coex_B = solve(E3_Pred_New2_B,H) % (c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) + O*b*c - 2*O*b*h*u)/(2*O*u*(c - h*u)) **OR** -(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) - O*b*c + 2*O*b*h*u)/(2*O*u*(c - h*u))
% The two possible sol'ns for each of A and B are the same, so still two
% solutions, and also means that it doesn't matter which equilibrium sol'n
% you use 
E3_Prey_New2_A = subs(E3_Prey,H,(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) + O*b*c - 2*O*b*h*u)/(2*O*u*(c - h*u)))
P3Equil_Coex_A = solve(E3_Prey_New2_A,P)
E3_Prey_New2_B = subs(E3_Prey,H,-(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) - O*b*c + 2*O*b*h*u)/(2*O*u*(c - h*u)))
P3Equil_Coex_B = solve(E3_Prey_New2_B,P)


% Jacobian matricies & dominant eigenvalues
d3PreydP = H*((2*O*P*c)/(O*h*P^2 + 1) - (2*O^2*P^3*c*h)/(O*P^2*h + 1)^2)
d3PreydH = (O*P^2*c)/(O*h*P^2 + 1) - u
d3PreddP = r - (2*H*O*P)/(O*h*P^2 + 1) + (2*H*O^2*P^3*h)/(O*P^2*h + 1)^2
d3PreddH = -(O*P^2)/(O*h*P^2 + 1)

% Invasion
J3Inv = subs([d3PreydP,d3PreydH,d3PreddP,d3PreddH],[P,H],[0,b/u])
% [ 0, -u, r, 0]
J3InvEig = r

% Coexistance
J3Coex_A = subs([d3PreydP,d3PreydH,d3PreddP,d3PreddH],[P,H],[0,(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) + O*b*c - 2*O*b*h*u)/(2*O*u*(c - h*u))])

J3Coex_B = subs([d3PreydP,d3PreydH,d3PreddP,d3PreddH],[P,H],[0,-(c*(O*(O*b^2 - 4*h*r^2*u^2 + 4*c*r^2*u))^(1/2) - O*b*c + 2*O*b*h*u)/(2*O*u*(c - h*u))])

% This is just becoming a biologically uninterpretable mess. We'll have to
% do this numerically.
