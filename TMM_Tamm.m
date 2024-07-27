%   Physical parameters
%       Central frequency wavelength
lambda0 = 300;  %   unit: um

%       Refraction index
n_Au = 5.4-1.75i;
n_Si = 3.4;
n_Vac = 1.;

%       Thickness
d_Au = 0.1; %   unit: um
d_Si = 0.75*lambda0/n_Si;
d_Vac = 0.25*lambda0/n_Vac;

%       DBR mirror multilayer total number
DBR_layer = 5;

%   Wavenumbers
k = 2*pi*linspace(0.6/(lambda0),1.4/(lambda0),1001);

%   Calculate DBR and Tamm cavity reflection coefficients for different
%   wavenumbers
reflections_Tamm = [];
reflections_DBR = [];

for i=1:length(k)
reflections_Tamm = [reflections_Tamm,Tamm_reflection(k(i),d_Si,n_Si,d_Vac,n_Vac,d_Au,n_Au,DBR_layer)];
reflections_DBR = [reflections_DBR,DBR_reflection(k(i),d_Si,n_Si,d_Vac,n_Vac,DBR_layer)];
end

%   Calculate Tamm mode field distribution
d = [d_Si,d_Vac,d_Si,d_Vac,d_Si,d_Au];
n = [n_Si,n_Vac,n_Si,n_Vac,n_Si,n_Au];
field = multilayer_field_distribution(k(500),d,n,100);

%   Plot reflection coefficients vs. frequency
figure
hold on

plot(k*lambda0/2/pi,reflections_Tamm)
plot(k*lambda0/2/pi,reflections_DBR)

legend('Tamm reflection','DBR reflection')

xlabel('Frequency(THz)')
ylabel('Reflection')

hold off

%   Plot Tamm mode field frequency
figure
hold on
x = field(1,:);
E = field(2,:);
plot(x,E)

%legend('Tamm reflection','DBR reflection')

%xlabel('Frequency(THz)')
%ylabel('Reflection')

hold off

    %   Functions
%   Transfer matrix formula (Macleod, H. A. (Hugh A. (2001). Thin-film optical filters / H.A. Macleod. (Third edition.). Institute of Physics Pub.)
function M = transfer_matrix(k,d,n)
    delta = k*n*d;
    M = [cos(delta),sin(delta)*n*1i;sin(delta)/n*1i,cos(delta)];
end

%   Reflection coefficient for multilayer structure
function r = multilayer_reflection(k,d,n)

    %   Error message when giving different sizes of thicknesses d and
    %   refractive index n.
    N = length(d);
    if length(d)~=length(n)
        prompt={'Mismatch size of n and d'};
        return
    end
    
    %   Calculate transfer matrix M, iterating for different layers.
    E_H_exit = [1;1]; % E and H fields at the end of the multilayer
    M = [1,0;0,1];
    for m=1:N
        M = M * transfer_matrix(k,d(m),n(m));
    end
    
    %   Calculate field amplitude reflection coefficient
    E_H_in = M * E_H_exit;
    C_B = E_H_in(2) / E_H_in(1);
    r = (1-C_B) / (1+C_B);
end

%   Electric field amplitude spatial distribution for multilayer structure
function field = multilayer_field_distribution(k,d,n,N_mesh)

    %   Error message when giving different sizes of thicknesses d and
    %   refractive index n.
    N = length(d);
    if length(d)~=length(n)
        prompt={'Mismatch size of n and d'};
        return
    end
    
    %   Calculate transfer matrix M, iterating for different layers.
    E_H_exit = [1;1]; % E and H fields at the end of the multilayer
    x = [0];
    x_present = 0;
    E = [1];
    E_H = E_H_exit;

    for m=1:N
        for j=1:N_mesh
            M = transfer_matrix(k,d(m)/N_mesh,n(m));
            E_H = M * E_H;
            x_present = x_present + d(m)/N_mesh;
            x = [x,x_present];
            E = [E, abs(E_H(1))];
        end
    end
    field = [x;E];
end

function R = DBR_reflection(k,d1,n1,d2,n2,N)
    
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

    %   Calculate power reflection of DBR mirror
    r = multilayer_reflection(k,d,n);
    R = r*conj(r);
end

function R = Tamm_reflection(k,d1,n1,d2,n2,d_metal,n_metal,N)
    %   Define list of thickness d and refractive index n
    d = [];
    n = [];

    %   DBR parameter
    for m=1:N
        if mod(m,2)==1
            d = [d,d1];
            n = [n,n1];
        else
            d = [d,d2];
            n = [n,n2];
        end
    end

    %   Metal layer parameter
    d = [d,d_metal];
    n = [n,n_metal];
    
    %   Power reflection
    r = multilayer_reflection(k,d,n);
    R = r*conj(r);
end