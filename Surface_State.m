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
DBR = make_DBR_unit_cells(d_Vac,n_Vac,d_Si,n_Si);

w = 2*pi* linspace(0.,1.8,100)/lambda0;
K = w;
t = t_analysis(K,DBR);

a = lambda0/4/n_Vac + lambda0*3/4/n_Si;
k = bandstructure(w,a,t);
plot(k/2/pi*a,w/2/pi*lambda0);
xlabel('k (2\pi/a)')
ylabel('Frequency(THz)')



function DBR = make_DBR_unit_cells(d1,n1,d2,n2)
    
    %   Define list of thickness d and refractive index n
    d = [];
    n = [];

    %   Assign element value, first layer is material 1 (d1,n1), second
    %   layer is material 2 (d2,n2), and so on.
    for m=1:3
        if mod(m,2)==1
            d = [d,d1/2.];
            n = [n,n1];
        else
            d = [d,d2];
            n = [n,n2];
        end
    end
    DBR = [d;n];
end

%   Transmission coefficient calculation
function t = t_analysis(K,parameters)
    
    %   Giving thickness and refractive index
    d = parameters(1,:);
    n = parameters(2,:);
    
    t = [];
    for j = 1:length(K)
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
            M = transfer_matrix(K(j),d(ms),n(ms));
            E_H = M * E_H;
        end
        field = [x;E/E(length(E))];
        
        %   Calculating reflection coefficients
        n_eff = E_H(2) / E_H(1);
        r = (1-n_eff) / (1+n_eff);
        R = abs(r)^2;
        T = 1-R;
        t = [t,(T)^0.5];
    end
end

%   Transfer matrix formula (Macleod, H. A. (Hugh A. (2001). Thin-film optical filters / H.A. Macleod. (Third edition.). Institute of Physics Pub.)
function M = transfer_matrix(k,d,n)
    delta = k*n*d;
    M = [cos(delta),sin(delta)/n*1i;sin(delta)*n*1i,cos(delta)];
end

function k = bandstructure(w,a,t)
    k = [];
    for i = 1:length(w)
        K = w(i);
        c_k = cos(K*a)/t(i);
    
        %ex ^ 2 - 2 * c_k * ex + 1 = 0;
        ex = (2 * c_k + sqrt(4 * c_k ^2 - 4))/2;
        k = [k, log(ex)/1i/a];
    end
end
