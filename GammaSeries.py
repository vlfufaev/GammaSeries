from sympy import *
import math

# AUXILIARY VALUES
z, p, Q, P, u, expr, expr1, U, Eq = symbols('z p Q P u expr expr1 U Eq ');

#PARAMETERS OF PROBLEM
N = '3';
M = '3';
q = symbols('q:' + M, real = True);
Q = 0;
for i in range(0, int(M)):
	Q = Q+ q[i]*(z+1)**(i + 1)/math.factorial(i+1);

c = symbols('c:' + N, real = True);
P = -1 + c[0]*p;
for i in range(1, int(N)):
	P = P - c[i]*p**(i + 1)*I;
	
#TAYLOR SERIES AFTER SUBSTITUTION t = 1 - u**2
expr = ((Q.subs(z, P.subs(p, p*(1 - u**2))) - Q.subs(z, P) -  (-1)*p*c[0]*q[0]*u**2)/( (-1)*p*c[0]*q[0]*u**2));
expr1 = 1;
for i in range(1, int(N)):
        expr1 = expr1 + ((expr)**i)*((-1)**i)*math.factorial(2*i)/(1 - 2*i)/(math.factorial(i))**2/(4**i);
        
#INTEGRAND
U = I*(u**2)*(expr1*(diff(P, p)).subs(p, p*(1 - u**2))/c[0]);

#SYSTEM
Eq = integrate(U, (u, 0, 1)).series(p, 0, int(N));
for i in range(1, int(N)):
        print(solve(re((diff(Eq, p, i)).subs([(p, 0), (c[0], exp(-I*pi/6))])), c[i]))
        
del(q, c, z, p, Q, P, u, expr, expr1, U, Eq, N, M)
