el=-75e-3;
vth=-50e-3;
vmax=10e-3;
vreset=-80e-3;
deltath=2e-3;
gl=10e-9;
cm=100e-12;
a=2e-9;
b=0.02e-9;
tsra=200e-3;

dt=0.5e-4;
tmax=1.5;
tvec=0:dt:tmax;

vm = zeros(size(tvec));
isra=zeros(size(tvec));
current=zeros(size(tvec));
vm(1)=el;
isra(1)=0;
iapp=500e-12;
pulseIdxs = tvec >= 0.5 & tvec <= 1;
current(pulseIdxs) = iapp;

for i=1:length(tvec)-1

    %forward Euler method
    if vm(i)>vth %if there is a spike
       vm(i)=vreset; %reset the voltage
       isra(i)=isra(i)+b; 
    end

    fv= (gl*( el-vm(i)+deltath*exp((vm(i)-vth)/deltath) )-isra(i)+current(i)) /cm;
    dv=dt*fv;
    vm(i+1)=vm(i)+dv;
    isra(i+1)=isra(i)+dt*(a*(vm(i)-el)-isra(i))/tsra;

end

figure;
tl = tiledlayout(3,1);

nexttile
plot(tvec,vm*1e3)
ylabel("Voltage (mV)")

nexttile
plot(tvec,current*1e9)
ylabel("Applied Current (nA)")

nexttile
plot(tvec, isra*1e9)
ylabel("SRA Current (nA)")