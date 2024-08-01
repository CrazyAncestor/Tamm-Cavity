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
DBR_layer = 111;

%   Give frequencies and parameters.
w = 2*pi* linspace(0.,2,1000)/lambda0/n_Vac;
a = lambda0/4/n_Vac + lambda0*3/4/n_Si;

%   Calculate the band structure
k = bandstructure(w,a,d_Vac,n_Vac,d_Si,n_Si);
att = bandstructure(2*pi* 1/lambda0/n_Vac,a,d_Vac,n_Vac,d_Si,n_Si)

plot(k/2/pi*a,w/2/pi*lambda0);
xlabel('k (2\pi/a)')
ylabel('Frequency(THz)')

%   Calculate DBR field
DBR = make_DBR(d_Si,n_Si,d_Vac,n_Vac,DBR_layer);
reflections_DBR = [];

for i=1:length(k)
    [r_DBR, field_DBR] = TMM_analysis(w(i),DBR,0);
    [r_Tamm, field_Tamm] = TMM_analysis(w(i),Tamm,0);
    reflections_DBR = [reflections_DBR , abs(r_DBR)^2];
end

%   Calculate DBR field
[min_value, w_idx] = max(reflections_DBR);
[r_DBR, field_DBR] = TMM_analysis(w(w_idx),DBR,100);
attenuation = Get_Attenuation(w(k_idx),DBR,100)

%   Plot DBR field enhancement
plot_field('DBR field enhancement',field_DBR,DBR_layer,d_Si,d_Vac)


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
function [r,t] = rt_analysis(K,parameters)
    
    %   Giving thickness and refractive index
    d = parameters(1,:);
    n = parameters(2,:);
    
    t = [];
    r = [];
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
        r0 = (1-n_eff) / (1+n_eff);
        t0 = 2 / (E_H(2)+E_H(1));
        r = [r,r0];
        t = [t,t0];
    end
end

%   Transfer matrix formula (Macleod, H. A. (Hugh A. (2001). Thin-film optical filters / H.A. Macleod. (Third edition.). Institute of Physics Pub.)
function M = transfer_matrix(k,d,n)
    delta = k*n*d;
    M = [cos(delta),sin(delta)/n*1i;sin(delta)*n*1i,cos(delta)];
end

function k = bandstructure(w,a,d1,n1,d2,n2)
    DBR_unit_cell = make_DBR_unit_cells(d1,n1,d2,n2);
    [r,t] = rt_analysis(w,DBR_unit_cell);
    k = [];
    for i = 1:length(w)
        K = w(i);
        delta = angle(t(i));

        c_k = cos(K*a - delta)/abs(t(i));
        k = [k, acos(c_k)/a];
    end
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
    
    %   Calculating reflection coefficients
    n_eff = E_H(2) / E_H(1);
    r = (1-n_eff) / (1+n_eff);
    
    %   Giving solutions of field
    field = [flip(x);flip(E)/E(length(E))];
end

%   Get attenuation factor
function attenuation = Get_Attenuation(k,parameters,N_mesh)
    
    %   Giving thickness and refractive index
    d = parameters(1,:);
    n = parameters(2,:);
    
    %   Defining variables
    E_H_exit = [1;1]; % E and H fields at the end of the multilayer
    x = [];
    x_present = sum(d);
    E = [];
    E_H = E_H_exit;
    
    %   Calculating the E,H fields with transfer matrices
    N = length(d);
    for m=1:N
        ms = N-m+1;
        x0 = [];
        E0 = [];

        
        for j=1:N_mesh
            M = transfer_matrix(k,d(ms)/N_mesh,n(ms));
            E_H = M * E_H;
            x_present = x_present - d(ms)/N_mesh;
            x0 = [x0, x_present];
            E0 = [E0, abs(E_H(1))];
        end
        if mod(ms,2)==1
           [max_E,idx] = max(E0);
           E = [E,E0(idx)];
           x = [x,x0(idx)];
        end  
    end
    
    %   Giving solutions of field
    field = [flip(x);log(flip(E))];
    P = polyfit(flip(x),log(flip(E)),1);
    slope = P(1);
    attenuation = slope;
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

function plot_field(name,field,layer_num,d_Si,d_Vac)
    figure
    title(name)
    hold on
    x = field(1,:);
    E = field(2,:);
    block_x = [0];
    x_present = 0;
    Height = max(E);
    for i=1:layer_num
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
end
