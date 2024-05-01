%% goy shell model

clear all
close all

%% declarations
global k N delta img forc lambda

PowMax = 10;     % max order of structure function
N = 20;          % number of shells
lambda = 2.0;    % the intershell ratio, usually set to 2.0
k0 = 0.05;       % the wavenumber at n = 0
delta = 0.5;     % GOY parameter in the nonlinear interaction term
nu = 1.0e-05;    % viscosity
DT = 0.001;      % time step
Nstat = 10^7;    % number of time steps, should be high enough for good statistics     
update_dt = 10^5;% # time steps for updating dt (diagnose in tau.dat)
decay_mode = 0;  % if 1, stop forcing after given # of time steps

img = complex(0.,1.);           % imaginary unity   
forc = complex(1.e-02,1.e-02);  % forcing 

Stat = 0;
dt = DT;

%% Compute wavenumers
   k = k0*lambda.^(0:N-1); 

%% Integrating factors for the viscous term
   for j=1:N
      f1(j) = complex(exp(-0.5 * nu * dt *k(j)*k(j)));
      f2(j) = (1 - f1(j)) / (nu * k(j) * k(j));
   end

%% set the average velocity, mu, to zero
   mu(1:N) = 0.0;

%% initialize the velocity according to Kolmogorov scaling
   stat = 0;
   dtmp = (k/k(1)).^(-1/3);
   u1  = complex(dtmp,0.0);
   du1 = goy_rhs(u1, stat);
   % working velocity
   u = u1;
   
%% initialize averaged structure function and energies    
   Sp(1:PowMax,1:N) = 0; eps_av = 0 ; en_av = 0 ; eps_forc = 0;
  
%% you should now open files for writing
       
% init time
    t = 0.0;
   
%% loop Nstat times

for stat = 0:Nstat-1
    
  %% Update average of u_n in order to estimate the smallest time scale
  % In the integration, the time step dt must be chosen according to the
  % smallest time scale in the cascade
    mu = mu + abs(u);
  % update time  
    t = t + dt;
       
  % check each update_dt time steps whether dt has to be modified
    if( mod(stat,update_dt)==0 && (stat>0) )
       dtmp = 1./(k.*mu/update_dt);
       min_dtmp = min(dtmp);

       fout3 = fopen('tau.dat','w');
       fprintf(fout3,'%d %e %e \n',[0:N-1; k; dtmp]);      
       fclose(fout3);

       % modify t, using the smallest time scale in the game...
       dt = min_dtmp / 50.;
       fprintf('doing %d from %d, dt now = %e \n',stat, Nstat, dt);
       % adapt integrating factors   
       for j=1:N
          f1(j) = complex(exp(-0.5 * nu * dt *k(j)*k(j)));
          f2(j) = (1 - f1(j)) / (nu *k(j) *k(j));
       end
       % reset average velocity to zero
       mu(1:N) = 0.0;
     end
    
   %% update mean quantities: structure function, energy and dissipations
   %  dump the instantaneous values of energy and dissipation rate every so 
   %  many time steps (say 1024) on a file
   %  dump their averages every so many tme steps (say 4096) on a new file
   
   %% How to compute
   % the structure function 
   %   for j = 1:N
   %      Sp(:,j) = Sp(:,j) + abs(u(j)).^(0:(PowMax-1))';
   %   end
   % the energy
   %   en = sum( real(u .* conj(u)) );
   % the dissipation rate
   %   eps = nu*sum(k .*k .*  real(u .* conj(u)));  
   % the energy forced into shell 1  
   %   eps_forc = forc * conj(u(1)));
   % a snapshot of the velocty field     
   %   xx = linspace(0, 2*pi/k(1), 128);                
   %   for ix=1:128
   %      uu(ix) = 0;
   %      for j=1:N
   %         uu(ix) = uu(ix) + exp(img*k(j)*xx(ix))*u(j);
   %      end
   %   end
   % finally take the real part of uu
 
        
   
   
   %% Integrate ODEs using the scheme of Pisarenko et al., POF 5, 2533 (1983)    
    du2 = goy_rhs(u1, stat);
    u2  = f1.*u1 + f2.* (1.5.*du2 - 0.5.*du1);
    u1  = u2;
    du1 = du2;
    
    % working u
    u = u2;
    
  %% if in decay mode, switch off forcing (here after 10^5 time steps)
    if((stat > 10^5) && decay_mode) 
        if(forc > 0) fprintf('\n Stopped forcing at t:%e\n', t);end;
        forc = 0;
    end;
end
    


%% RHS of GOY model
function du = goy_rhs(u, stat)
global k N delta img forc lambda

% Evaluate r.h.s. of eqn. of motion

coeffb = -delta / lambda;
coeffc = -(1.0 - delta) / (lambda^2);

cu = conj(u);

du(1) = img*k(1)*(cu(2)*cu(3));
du(2) = img*k(2)*(cu(3)*cu(4)+coeffb*cu(1)*cu(3));

aa = cu(4:N-1).*cu(5:N);
bb = cu(2:N-3).*cu(4:N-1);
cc = cu(1:N-4).*cu(2:N-3);
du(3:N-2) = img.*k(3:N-2).*(aa+coeffb*bb+coeffc*cc);

du(N-1) = img*k(N-1)*(coeffb*cu(N-2)*cu(N) + coeffc*cu(N-3)*cu(N-2));
du(N)   = img*k(N)*(coeffc*cu(N-2)*cu(N-1));

% Force lowest shell
du(1) = du(1) + forc;


end
