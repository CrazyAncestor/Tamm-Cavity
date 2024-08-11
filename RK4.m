global e
global m
global hbar
global Tesla
global V_m
global Angstrom
global light_speed
global G_GaAs

e = 1;
m = 0.07;
hbar = 1/2/pi;
Tesla = 1/35.75;
V_m = 1/96.5e6;
Angstrom = 0.1/270;
light_speed = 3e6 * Angstrom;
G_GaAs = 2*pi/(5.7* Angstrom);

E0 = 0 * Tesla * light_speed/3.4;
B0 = 5 * Tesla;
Em = 2e0 * Tesla * light_speed/3.4;
Bm = 2e0 * Tesla;
omega = 10*2*pi; % Unit THz

k0 = [0,G_GaAs*5e-5];

fc = e*B0/m/2/pi
CR_motion([E0,B0,Em,Bm,omega],k0,50/fc,1e-3)

function CR_motion(field,k0,tf,dt)
    global E0
    global B0
    global Em
    global Bm
    global omega
    global Angstrom
    global G_GaAs

    E0 = field(1);B0 = field(2);Em = field(3);Bm = field(4); omega = field(5);
    N = int32(tf/dt);
    [t,k] = RK4_solution(@EoM_lowV,[0,[k0(1),k0(2)]],dt,N);
    r = motion2(t,k);
    plotxy(k,'Momentum',G_GaAs)
    plotxy(r,'Spatial Orbit',Angstrom)
    
    %   FFT of momentum
    spectrum_kx = abs(fft(k(1,:)));
    fs = 1/dt;
    N = length(spectrum_kx);
    f = (0:N-1)*fs/N;
    figure
    plot(f(1:int32(N*0.001)),spectrum_kx(1:int32(N*0.001)))

    %   FFT of momentum
    phase_x = [r(1,:);k(1,:)];
    plotxy(phase_x,'Phase Diagram x',1)
    phase_y = [r(2,:);k(2,:)];
    plotxy(phase_y,'Phase Diagram y',1)
end

function plotxy(xy,title_name,unit_length)
    figure
    title(title_name)
    hold on
    x = xy(1,:);y = xy(2,:);
    plot(x/unit_length,y/unit_length)
    grid("on")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    hold off
end

function r = motion(t,k)

    global e
    global m
    global hbar
    
    EB = EB_field(t);
    E = EB(1,:); B = EB(2,:);
    
    kx = k(1,:);ky = k(2,:);
    rx =  ky./B*(e*hbar);
    ry = -kx./B*(e*hbar) - E.*t./B;
    r = [rx;ry];
end

function r = motion2(t,k)

    global m
    global hbar
    
    kx = k(1,:);ky = k(2,:);
    rx = []; ry = [];
    x=0;y=0;
    dt = t(2)-t(1);

    for i = 1:length(t)
        x = x + hbar*kx(i)/m*dt;
        y = y + hbar*ky(i)/m*dt;
        rx = [rx,x];
        ry = [ry,y];
    end
    r = [rx;ry];
end

function f = test(x,y)

    f = [y(1) ; y(2)];
end

function EB = EB_field(t)
    global E0
    global B0
    global Em
    global Bm
    global omega
    
    N = length(t);
    E = E0*ones(1,N) + Em*sin(omega*t)+ Em*sin(0.318*t*2*pi);
    B = B0*ones(1,N) + Bm*sin(omega*t)+ Bm*sin(0.318*t*2*pi);

    EB = [E;B];
end

function dk_dt = EoM_lowV(t,k)

    global e
    global m
    global hbar

    EB = EB_field(t);
    E = EB(1);B = EB(2);
    
    kx = k(1);ky = k(2);
    dky_dt = -(e*B/m)*kx;
    dkx_dt =  (e*B/m)*ky + e*E/hbar;
    dk_dt = [dkx_dt;dky_dt];
end

function [xs,ys] = RK4_solution(func,boundary,dx,N)
    x0 = boundary(1); y0 = transpose(boundary(2:length(boundary)));
    x = x0;
    y = y0;
    xs = [x0];
    ys = [y0];
    
    for i=1:N
        k1 = func(x,y);
        k2 = func(x + dx/2.,y + dx/2.*k1);
        k3 = func(x + dx/2.,y + dx/2.*k2);
        k4 = func(x + dx,y + dx*k3);
        y = y + dx/6.*(k1 + k2 * 2 + k3 * 2 + k4);
        x = x + dx;
        xs = [xs,x];
        ys = [ys,y];
    end
end