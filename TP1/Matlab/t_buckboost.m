clear all;
close all;

rl = sym('rl');
rc = sym('rc');
C  = sym('C');
L  = sym('L');
D  = sym('D');
R  = sym('R');
vd = sym('vd');

A1 = [[-rl/L 0]; [0 -1/(C*(R+rc))]];
B1 = [1/L ; 0];
C1 = [0 R/(R+rc)];

% OFF State Same as Buck
A2 = 1/(R+rc)*[[-(R*rc+rl*rc+rl*R)/L -R/L]; [R/C -1/C]];
B2 = [0; 0];
C2 = 1/(R+rc)*[R*rc R];

A0 = A1*D+A2*(1-D);
B0 = B1*D+B2*(1-D);
C0 = C1*D+C2*(1-D);

disp('v0/vd steady state');
X_steady    = simplify(-inv(A0)*B0*vd);
t_vd_steady = simplify(-C0*inv(A0)*B0);
pretty(t_vd_steady)

rl = 0;
rc = 0;

disp('v0/vd steady state rl=0 rc=0');
t_vd_steady_r0 = simplify(subs(t_vd_steady));
pretty(t_vd_steady_r0)

clear rl, rc;
rl = sym('rl');
rc = sym('rc');

s = sym('s');
I = [[1 0];[0 1]];

t_vd = C0*inv(s*I-A0)*B0;

rc = 10e-3;
rl = 20e-3;
C  = 2000e-6;
L  = 5e-6;
R  = 200e-3;
D  = 0.625;
vd = 8;
fs = 200e3;

rl=0;
rc=0;

t_vd_n = simplify(subs(t_vd));
[n_vd, d_vd] = numden(t_vd_n);
n_vd_n = sym2poly(n_vd);
d_vd_n = sym2poly(d_vd);

% Normalization
n_vd_n = n_vd_n/d_vd_n(3);
d_vd_n = d_vd_n/d_vd_n(3);

disp('vo/vd');
sys_vd = tf(n_vd_n, d_vd_n)
figure(1);
set(1, 'NumberTitle', 'off'); 
set(1, 'Name', 'v0/vd Bode'); 

bode(sys_vd);
figure(2);
set(2, 'NumberTitle', 'off'); 
set(2, 'Name', 'v0/vd pzmap'); 
pzmap(sys_vd);

subs(X_steady)

t_d = C0*inv(s*I-A0)*((A1-A2)*X_steady + (B1-B2)*vd)+(C1-C2)*X_steady;
t_d_n = simplify(subs(t_d));

[n_d, d_d] = numden(t_d_n);
n_d_n = sym2poly(n_d);
d_d_n = sym2poly(d_d);

% Normalization
n_d_n = n_d_n/d_d_n(3);
d_d_n = d_d_n/d_d_n(3);

disp('vo/d');
sys_d = tf(n_d_n, d_d_n)
figure(3);
set(3, 'NumberTitle', 'off'); 
set(3, 'Name', 'v0/d Bode'); 

bode(sys_d);
figure(4);
set(4, 'NumberTitle', 'off'); 
set(4, 'Name', 'v0/d pzmap'); 
pzmap(sys_d);