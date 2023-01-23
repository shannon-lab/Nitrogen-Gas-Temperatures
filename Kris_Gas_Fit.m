%%%This code simulates the rovibrational bands of the second positive N2
%%%system, the v'−v" (0−0) transition
%%%Output is rotational temperatures with error bar
%%% For an overview of the physics contained within this script please read
%%% the thesis of Kristopher Ford, which is included in the repository.

clear all
%Raw data to be read, 1st column is wavelength (nm), 2nd set of 
%columns is intensity, 3rd set of columns is reference, and last column is 
%std_error. It is required that the data have a variable called 'alldata' which contains
%that format. an example would be 1000x10 for 1000 wavelength data points
mat_title = 'testdata.mat';
load(mat_title);

num_bootstrap=10; %number of times you want to re−simulate the T_r to get simulation error

% labeling the input data
[row col] = size(alldata);
expt_lam = alldata(:,1);
expt_int = alldata(:,2:(col-2)/2+1);
expt_ref = alldata(:,(col-2)/2+2:col-1);
sigma = alldata(:,col);

%%
% find the indices for truncating the data to the (0,0) v',v" transition
left = find(abs(expt_lam-335.7)<0.008);
right = find(abs(expt_lam-337.3)<0.008);
expt_lam = expt_lam(left:right,:);
expt_int = expt_int(left:right,:);
expt_ref = expt_ref(left:right,:);%added by Josh M.
sigma = sigma(left:right,:);
rel_expt_err = sigma./(mean(expt_int,2)-mean(expt_ref,2));
%%
%Total number of rotational lines for each PQR branch; could change J's but recommend 40.
Jmax = 40;
%Rotational temperatures to attempt. The script will change the T−rotation limits if convergence isn't found
%Users can change this to a starting range that suits the experiment, wide
%ranges with small step size take longer to run.
T_rotation = linspace(350,550,ceil((550-350)/2)+1);% in K with a spacing of 1K per step

%%% SCALING FACTORS %%%
%change this value to match the FWHM of your spectrometer
fwhm = 0.031; %fwhm on czerny turner spectrometer, nm

gauss_sig = fwhm/(2*(2*log(2))^0.5); %Do not change

%%% CONSTANTS %%%
k_b = 1.38E-23; % Boltzman constant (J/K)
c = 3.0*10^8; % speed of light (m/s)
h = 6.626070*10^-34; % planks constant in (J*s)
q_franck = 0.4515; %Franck−Condon factor


%Rotational C3pi_u constants are in 1st column, B3pi_g in second
%These are from Roux 1988, https://doi.org/10.1139/p89−023
Bv = [1.8153 1.62872]; %(cm−1)
A = [39.134 42.234]; %(cm−1)
Y = A./Bv;
Dv = [0.595E-5 0.581E-5]; %(cm−1) called 'D' in Roux paper but Dv in textbooks
waveorigin = 29670.942; %the (0−0) transition origin from Roux 1988 (cm−1)

%Rotational wavenumber calculations from Budo 1933, reproduced in
%Herzberg's Spectra of Diatomic Molecules Vol.1 pg 235
for J=0:Jmax
for i=1:2
    jay = J+1;
    Z1 = Y(i)*(Y(i)-4)+4/3+4*J*(J+1);
    Z2 = 1/3/Z1*(Y(i)*(Y(i)-1)-4/9-2*J*(J+1));
    F1(jay,i)=Bv(i)*(J*(J+1)-Z1^0.5-2*Z2)-Dv(i)*(J-1/2)^4;
    F2(jay,i)=Bv(i)*(J*(J+1)+4*Z2)-Dv(i)*(J+1/2)^4;
    F3(jay,i)=Bv(i)*(J*(J+1)+Z1^0.5-2*Z2)-Dv(i)*(J+3/2)^4;
    %F1,2,3 are now wavenumbers for the B state in column2, C in column1
    %F(1,i) is F(J=0,i) and F(Jmax+1,i) is F(J=40)
end
end

jaymax=Jmax+1;
%Calculating difference in term values for the wavelength origins
wave_p1 = F1(3:jaymax-1,1)-F1(4:jaymax,2);
wave_p2 = F2(3:jaymax-1,1)-F2(4:jaymax,2);
wave_p3 = F3(3:jaymax-1,1)-F3(4:jaymax,2);
waveP = [wave_p1;wave_p2;wave_p3]+waveorigin; %waveorigin is the bandhead
lamP = waveP.^-1*10^7; %conversion to nm

wave_q1 = F2(3:jaymax,1)-F2(3:jaymax,2);
wave_q2 = F3(2:jaymax,1)-F3(2:jaymax,2);
waveQ = [wave_q1;wave_q2]+waveorigin;
lamQ = waveQ.^-1*10^7;

wave_r1 = F1(2:jaymax,1)-F1(1:jaymax-1,2);
wave_r2 = F2(3:jaymax,1)-F2(2:jaymax-1,2);
wave_r3 = F3(4:jaymax,1)-F3(3:jaymax-1,2);
waveR = [wave_r1;wave_r2;wave_r3]+waveorigin;
lamR = waveR.^-1*10^7;
%Combining the P,Q,R branches into one simulated lambda vector
lambda=[lamP;lamQ;lamR];

%Line Strength factors from Kovacs 1969 pg 132
%u and C expressions
Jvec = linspace(0,jaymax,Jmax+2)';
why = Y; %Making Y single valued so this is easier to debug
Y=Y(1);

%Calculating the C state constants
u1(:,1) = (Y*(Y-4)+4*Jvec.^2).^0.5;
u3(:,1) = (Y*(Y-4)+4*(Jvec+1).^2).^0.5;
C1(:,1) = Y*(Y-4)*Jvec.*(Jvec+1) + 2*(2*Jvec+1).*(Jvec-1).*Jvec.*(Jvec+1);
C2(:,1) = Y*(Y-4) + 4*Jvec.*(Jvec+1);
C3(:,1) = Y*(Y-4)*(Jvec-1).*(Jvec+2) + 2*(2*Jvec+1).*Jvec.*(Jvec+1).*(Jvec+2);
u1_plus(:,1) = u1(:,1)+(Y-2);
u1_minus(:,1) = u1(:,1)-(Y-2);
u3_plus(:,1) = u3(:,1)+(Y-2);
u3_minus(:,1) = u3(:,1)-(Y-2);

%Switching Y to Y(B state) and calculating those constants
Y = why(2);
u1(:,2) = (Y*(Y-4)+4*Jvec.^2).^0.5;
u3(:,2) = (Y*(Y-4)+4*(Jvec+1).^2).^0.5;
C1(:,2) = Y*(Y-4)*Jvec.*(Jvec+1) + 2*(2*Jvec+1).*(Jvec-1).*Jvec.*(Jvec+1);
C2(:,2) = Y*(Y-4) + 4*Jvec.*(Jvec+1);
C3(:,2) = Y*(Y-4)*(Jvec-1).*(Jvec+2) + 2*(2*Jvec+1).*Jvec.*(Jvec+1).*(Jvec+2);
u1_plus(:,2) = u1(:,1)+(Y-2);
u1_minus(:,2) = u1(:,1)-(Y-2);
u3_plus(:,2) = u3(:,1)+(Y-2);
u3_minus(:,2) = u3(:,1)-(Y-2);

%Kovacs expressions, easier to work from the 1966 Paper where Lambda=1
%rather than the textbook which gives formulation for any Lambda. Note that
%this is the intermediate hund's coupling expression for the main P,Q,R
%lines. ie. P_1, P_2, P_3. Don't use the P_Q_12 stuff
%Note that J refers to the upper level, J', and jay must be used as an index
%in order to get a J=0 calculation later.
Y=why;
for J = 2:Jmax-1
    jay = J+1;
    S_P1(jay) = (J-1)*(J+1)/(16*J*C1(jay-1,1)*C1(jay,2)) * ...
    (8*(J-2)*(J-1)*J*(J+1)+(J-2)*(J+2)*u1_minus(jay-1,1) * ...
    u1_minus(jay,2)+J^2*u1_plus(jay-1,1)*u1_plus(jay,2))^2;
    S_P2(jay) = (J-1)*(J+1)/(J*C2(jay-1,1)*C2(jay,2)) * ...
    ((Y(1)-2)*(Y(2)-2)+2*((J-2)*(J+2)+J^2))^2;
    S_P3(jay) = (J-1)*(J+1)/(16*J*C3(jay-1,1)*C3(jay,2)) * ...
    (8*(J-1)*J*(J+1)*(J+2)+(J-2)*(J+2)*u3_plus(jay-1,1) * ...
    u3_plus(jay,2)+J^2*u3_minus(jay-1,1)*u3_minus(jay,2))^2;
end

S_P1 = S_P1(3:jaymax-1);
S_P2 = S_P2(3:jaymax-1);
S_P3 = S_P3(3:jaymax-1);


for J=1:Jmax
    jay=J+1;
    S_Q1(jay) = (J-1)^2*(2*J+1)/(4*J*(J+1)*C1(jay,1)*C1(jay,2)) * ...
    (4*(J-1)*(J+1)^2+(J+2)*u1_minus(jay)*u1_minus(jay,2))^2;
    S_Q2(jay) = (2*J+1)/(J*(J+1)*C2(jay,1)*C2(jay,2)) * ...
    ((Y(1)-2)*(Y(2)-2)+2*((J-2)*(J+2)+J^2))^2;
end
S_Q1 = S_Q1(3:jaymax);
S_Q2 = S_Q2(2:jaymax);

for J=1:Jmax
    jay=J+1;
    S_R1(jay) = J*(J+2)/(16*(J+1)*C1(jay+1,1)*C1(jay,2)) * ...
    (8*(J-1)*J*(J+1)*(J+2)+(J-1)*(J+3)*u1_minus(jay+1,1) * ...
    u1_minus(jay,2)+(J+1)^2*u1_plus(jay+1,1)*u1_plus(jay,2))^2;
    S_R2(jay) = J*(J+2)/((J+1)*C2(jay+1,1)*C2(jay,2)) * ...
    ((Y(1)-2)*(Y(2)-2)+2*((J-1)*(J+3)+(J+1)^2))^2;
	S_R3(jay) = J*(J+2)/(16*(J+1)*C3(jay+1,1)*C3(jay,2)) * ...
    (8*J*(J+1)*(J+2)*(J+3)+(J-1)*(J+3)*u3_plus(jay+1,1) * ...
    u3_plus(jay,2)+(J+1)^2*u3_minus(jay+1,1)*u3_minus(jay,2))^2;

end

S_R1 = S_R1(2:jaymax);
S_R2 = S_R2(3:jaymax);
S_R3 = S_R3(4:jaymax);

SP = [S_P1 S_P2 S_P3]';
SQ = [S_Q1 S_Q2]';
SR = [S_R1 S_R2 S_R3]';
S_all = [S_P1 S_P2 S_P3 S_Q1 S_Q2 S_R1 S_R2 S_R3]';

timestep = 1;

%randomly select intensity and reference sets then generate new spectra
%called 'signal'
int_resample = randi([1 (col-2)/2],1,num_bootstrap);
ref_resample = randi([1 (col-2)/2],1,num_bootstrap);
for i=1:num_bootstrap
    signal(:,i)=expt_int(:,int_resample(i)) - expt_ref(:,ref_resample(i));
end

%%

for boot = 1:num_bootstrap
	boot_int = signal(:,boot);

    for iteration = 1:length(T_rotation)
        %%% TEMPERATURES %%%
        T_rot = T_rotation(iteration); % in K
	
        %multiplication factor for F in exp(F/k_B/T_r) expression. Note that THIS
        %is what is changing between iterations to fit T_r. 100 factor is to change
        %from cm^−1 to m^−1 and work with the h,c constants
        mult = -100*c*h/k_b/T_rot;

        expP1 = exp(F1(3:jaymax-1,1)*mult);
        expP2 = exp(F2(3:jaymax-1,1)*mult);
        expP3 = exp(F3(3:jaymax-1,1)*mult);
	
        expQ1 = exp(F2(3:jaymax,1)*mult);
        expQ2 = exp(F3(2:jaymax,1)*mult);
	
        expR1 = exp(F1(2:jaymax,1)*mult);
        expR2 = exp(F2(3:jaymax,1)*mult);
        expR3 = exp(F3(4:jaymax,1)*mult);
        exponent = [expP1; expP2; expP3; expQ1; expQ2; expR1; expR2; expR3];

        Int_stems = 1./lambda.^4*q_franck.*S_all.*exponent;
        dim = length(lambda); %size of the 'pure line' spectrum
        data_size = length(boot_int);
	
        %Fitting an oversampled simulation (higher resolution)
        Int_10xfinal = zeros(data_size*10,1);
        tenxlam = linspace(expt_lam(1),expt_lam(data_size),data_size*10)';
	
        %Summing together all the Gaussians for every single line
        for i = 1:data_size*10
            for j = 1:dim-1
                Int_10xfinal(i) = Int_10xfinal(i) + ...
                Int_stems(j)*normpdf(tenxlam(i),lambda(j),gauss_sig);
            end
	
        end
        Int_10xfinal = Int_10xfinal/max(Int_10xfinal)*max(boot_int);
	
	
        %shift the simulated data in the x−axis to get alignment at the
        %intensity maxima for simulation and experiment
        [discard, index] = max(Int_10xfinal);
        [discard, expt_index] = max(boot_int);
        shift= tenxlam(index)-expt_lam(expt_index); %shift in nm
        % shift = 0; %testing out no shift
        expt_lam_final = expt_lam+shift;
	
        %Interpolate and find error for simulation
        y_sim = interp1(tenxlam, Int_10xfinal, expt_lam_final);
        y_sim(isnan(y_sim))=0; %remove NaN results from interpolation function
        chi_error(iteration) = sum((((y_sim-boot_int)./sigma)).^2);
        %chi_error(iteration) = sum((((y_sim-boot_int))).^2); %switch to least squares
        relative_err(iteration) = sum(abs(y_sim-boot_int)./boot_int);
        allshift(iteration) = shift;
        all_expt_lam(:,iteration) = expt_lam_final;
	
        all_int_final(:,iteration) = y_sim;
	
    end
	
	%select best simulation and check if chi is minimized at the limit.
	%If so, change T_rotation limits and re-run simulations to find min
	[discard best_iteration] = min(chi_error);
	if best_iteration == 1
        T_rotation=T_rotation-150;
        boot=boot-1;
        fprintf('FAILINGDOWN')
	elseif best_iteration == length(T_rotation)
        T_rotation = T_rotation+150;
        boot=boot-1;
        fprintf('FAILINGUP')
    else
        bootstrap_temp(boot) = T_rotation(best_iteration);
        fprintf(['working' num2str(boot)])
    end
	
	if T_rotation(best_iteration) > 1000 || T_rotation(best_iteration) < 200
        temperature_limit_exceeded
    end
	
end
%%
%Average the bootstrap temperature results, find the simulation error
final_Tr = mean(bootstrap_temp(find(bootstrap_temp)));
fprintf('\n')
disp(['Final_Tr = ' num2str(final_Tr)])
rel_sim_err = std(bootstrap_temp(find(bootstrap_temp)))/final_Tr;
%overallerr will ignore any (inf) values that may inadvertently be a part
%of rel_expt_err in the sampled range
%rel_expt_err = rel_expt_err(left:right,1); %modifies error to the range simulated
overallerr = sqrt(mean(rel_expt_err(~isinf(rel_expt_err)))^2+rel_sim_err^2)*final_Tr;
fprintf('\n')
disp(['Final_err = ' num2str(overallerr)])
fprintf('\n')
disp(['Best chi^2 = ' num2str(min(chi_error))])
bestsimulation_int = all_int_final(:,best_iteration);

%Plot the chi−squared vs temperature for last run
figure();
plot(T_rotation,chi_error,'o')
title('Test Statistic vs. Simulation Temperature')
ylabel('\chi^2 (unitless)')
xlabel('T_r (K)')


%Plot the bootstrap results for gas temperature
figure();
plot(bootstrap_temp)
xlabel('Bootstrap Iteration Number')
ylabel('T_r Result')

% %Plot the relative error of the experimental data
% figure();
% plot(alldata(:,1),rel_expt_err,'ro')
% ylabel('Sigma/Avg_I (%)')
% xlabel('Wavelength (nm)')
% title('Experimental Error')
%%

%Plot the shifted experimental data and the simulation lines
figure()
plot(expt_lam_final,expt_int(:,1),'r',expt_lam_final,bestsimulation_int,'b')        
title(strcat('T_r=', num2str(final_Tr), ' FWHM=', num2str(fwhm)))
legend('Experiment','Simulation','Location','northwest')
ylabel('Intensity (a.u.)')
xlabel('Wavelength (nm)')
xlim([336 337.5])

%plot unshifted data
% figure()
% plot(expt_lam,expt_int,'r',sim_lam,Int_final,'b')
% title(strcat('T_r=', num2str(T_rot), ' FWHM=', num2str(fwhm)))
% legend('Experiment','Simulation')
% ylabel('Intensity (a.u.)')
% xlabel('Wavelength (nm)')
% xlim([335.5 337.5])

% % Plot all the simulation lines
% figure()
% stem(lambda,Int_stems,'.')
% xlim([335.5 337.5])
% title(strcat('Pure Simulation Lines Tr=', num2str(T_rot)))
% ylabel('Intensity (a.u.)')
% xlabel('Wavelength (nm)')        
