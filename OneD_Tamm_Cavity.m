%   Physical parameters
%       Central frequency wavelength
lambda0 = 300;  %   unit: um

%       Refraction index
n_Au = 400-500i;
n_Si = 3.4;
n_Vac = 1.;

%       Thickness
d_Au = 0.1; %   unit: um
d_Si = 0.75*lambda0/n_Si;
d_Vac = 0.25*lambda0/n_Vac;

%       DBR mirror multilayer total number
DBR_layer = 5;

%   Build DBR & Tamm cavity parameters
DBR = make_DBR(d_Si,n_Si,d_Vac,n_Vac,DBR_layer);
Tamm = make_Tamm(DBR, d_Au, n_Au);

%   Wavenumbers
k = 2*pi*linspace(0.6/(lambda0),1.4/(lambda0),10001);

%   Calculate DBR and Tamm cavity reflection coefficients for different
%   wavenumbers
reflections_Tamm = [];
theta_Tamm = [];
reflections_DBR = [];
theta_DBR = [];

for i=1:length(k)
    [r_DBR, field_DBR] = TMM_analysis(k(i),DBR,0);
    [r_Tamm, field_Tamm] = TMM_analysis(k(i),Tamm,0);
    reflections_DBR = [reflections_DBR , abs(r_DBR)^2];
    reflections_Tamm = [reflections_Tamm , abs(r_Tamm)^2];
    theta_Tamm = [theta_Tamm,angle(r_Tamm)];
    theta_DBR = [theta_DBR,angle(r_DBR)];
end

%   Calculate Tamm mode field enhancement
[min_value, k_idx] = min(reflections_Tamm);
[r_Tamm, field_Tamm] = TMM_analysis(k(k_idx),Tamm,100);

%   Plot reflection coefficients vs. frequency
figure
title('Reflection coefficient vs. frequency')
hold on

plot(k*lambda0/2/pi,reflections_Tamm)
plot(k*lambda0/2/pi,reflections_DBR)

legend('Tamm reflection','DBR reflection')

xlabel('Frequency(THz)')
ylabel('Reflection')

hold off

%   Smith Chart
Smith_Chart(reflections_DBR,theta_DBR,'Smith Chart (DBR cavity)')
Smith_Chart(reflections_Tamm,theta_Tamm,'Smith Chart (Tamm cavity)')

%   Plot Tamm mode field enhancement
figure
title('Tamm mode enhancement')
hold on
x = field_Tamm(1,:);
E = field_Tamm(2,:);
block_x = [0];
x_present = 0;
Height = max(E);
for i=1:DBR_layer
    if mod(i,2)==1
        block_x=[block_x,d_Si+x_present]; 
        x_present = x_present + d_Si;
        patch('Faces',[1 2 3 4],'Vertices',[ block_x(i) 0; block_x(i+1) 0; block_x(i+1) Height; block_x(i) Height],'FaceColor',[198 198 198]./255)
    else
        block_x=[block_x,d_Vac+x_present];
        x_present = x_present + d_Vac;
        patch('Faces',[1 2 3 4],'Vertices',[ block_x(i) 0; block_x(i+1) 0; block_x(i+1) Height; block_x(i) Height],'FaceColor',[255 255 255]./255)
    end
    
end
plot(x,E)
xlabel('position(um)')
ylabel('Electric Field Enhancement')

hold off

    %   Functions
%   Transfer matrix formula (Macleod, H. A. (Hugh A. (2001). Thin-film optical filters / H.A. Macleod. (Third edition.). Institute of Physics Pub.)
function M = transfer_matrix(k,d,n)
    delta = k*n*d;
    M = [cos(delta),sin(delta)/n*1i;sin(delta)*n*1i,cos(delta)];
end

%   Electric field amplitude spatial enhancement for multilayer structure
function [r, field] = TMM_analysis(k,parameters,N_mesh)
    
    %   Giving thickness and refractive index
    d = parameters(1,:);
    n = parameters(2,:);
    
    %   Defining variables
    E_H_exit = [1;1]; % E and H fields at the end of the multilayer
    x = [sum(d)];
    x_present = sum(d);
    E = [1];
    E_H = E_H_exit;
    
    %   Calculating the E,H fields with transfer matrices
    N = length(d);
    for m=1:N
        ms = N-m+1;
        if N_mesh==0
            M = transfer_matrix(k,d(ms),n(ms));
            E_H = M * E_H;
        else
            for j=1:N_mesh
                M = transfer_matrix(k,d(ms)/N_mesh,n(ms));
                E_H = M * E_H;
                x_present = x_present - d(ms)/N_mesh;
                x = [x, x_present];
                E = [E, abs(E_H(1))];
            end
        end
    end
    field = [x;E/E(length(E))];
    
    %   Calculating reflection coefficients
    n_eff = E_H(2) / E_H(1);
    r = (1-n_eff) / (1+n_eff);

end

function DBR = make_DBR(d1,n1,d2,n2,N)
    
    %   Define list of thickness d and refractive index n
    d = [];
    n = [];

    %   Assign element value, first layer is material 1 (d1,n1), second
    %   layer is material 2 (d2,n2), and so on.
    for m=1:N
        if mod(m,2)==1
            d = [d,d1];
            n = [n,n1];
        else
            d = [d,d2];
            n = [n,n2];
        end
    end
    DBR = [d;n];
end

function Tamm = make_Tamm(DBR,d_metal,n_metal)
    d = DBR(1,:);
    n = DBR(2,:);
    d = [d,d_metal];
    n = [n,n_metal];
    Tamm = [d;n];
end

function Smith_Chart(amp,ang,title_name)
    figure
    title(title_name)
    hold on
    x_smith = [];
    y_smith = [];
    for i=1:length(ang)
        x_smith = [x_smith,amp(i)^0.5*cos(ang(i))];
        y_smith = [y_smith,amp(i)^0.5*sin(ang(i))];
    end
    plot(x_smith,y_smith);
    hold off
    grid("on")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    axis([-1 1 -1 1])
end