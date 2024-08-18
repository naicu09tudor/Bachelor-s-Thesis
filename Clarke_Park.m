% obtaining three phase voltages
% transforming two phase, stationary and rotary 
fs=500;
Ts=1/fs;
t=0:Ts:1;
f=2;
w=2*pi*f;
Vm=220;
Va=Vm*cos(w*t);
Vb=Vm*cos(w*t-2*pi/3);
Vc=Vm*cos(w*t-4*pi/3);

% first alternative for abc to alpha-beta
% K=(2/3)*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2];
Val=2/3*(Va-0.5*Vb-0.5*Vc);
Vbet=2/3*(sqrt(3)/2*Vb-sqrt(3)/2*Vc);
% second alternative for abc to alpha-beta
V=(2/3)*(Va+Vb*exp(j*2*pi/3)+Vc*exp(-j*2*pi/3)); % voltage vector
Valpha=real(V);
Vbeta=imag(V);
% alpha-beta to dq
Vd=Valpha.*cos(w*t)+Vbeta.*sin(w*t);
Vq=-Valpha.*sin(w*t)+Vbeta.*cos(w*t);


%plot
subplot(2,2,1)
plot(t,Va,t,Vb,t,Vc);
grid
subplot(2,2,2)
plot(t,Valpha,t,Vbeta);
grid
subplot(2,2,3)
plot(t,Vd,t,Vq);
grid
subplot(2,2,4)
plot(Valpha/Vm,Vbeta/Vm)
grid
