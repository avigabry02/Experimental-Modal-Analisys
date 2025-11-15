clc
close all
clear all 


%% Loading of the files

% selecting main folder where are all the data subfolder sets
mainFolder = '/Users/gabri/Documents/università /PoliMi/Lessons/First Year/Second Semester/Data Analysis for Mechanical Systems Identification/Project/MatlabModes/test modal analysis';
% mainFolder = 'C:\Users\Utente\Desktop\LAUREA_MAGISTRALE\SEMESTER_2\DATA ANALYSIS FOR MECHANICAL SYSTEM IDENTIFICATION\PROJECT\experimental_data';
% mainFolder = 'C:\Users\Carlo Tundo\Desktop\PoliMi\Master 1st year\DATA ANALYSIS\assignement\modal analysis\test modal analysis';
folderList = genpath(mainFolder);           % Get all paths
folderListCell = strsplit(folderList, pathsep);  % Split into individual folders
folderListCell = folderListCell(~cellfun('isempty', folderListCell));  % Remove empty entries

folderListCell(1) = []; % deleting main folder path 

% definition of the dimension
ntest = 10;
nset = 5;
test = cell(nset, 1);

for ii = 1:nset % for each set 
    % preallocation 
    test{ii}.data = cell(ntest, 1);
    % selceting folder 
    currentFolder = folderListCell{ii};
    files = dir(fullfile(currentFolder, '*.mat'));
    
    % acquisition data
    for jj = 1:length(files)
        filePath = fullfile(currentFolder, files(jj).name);
        test{ii}.data{jj} = load(filePath);
    end
end

n_acc = size(test{1}.data{1}.Dati, 2) - 1;

%% Sensitivity of the dataset and creation of time vector
hammer_sens_steel = 2.488; % [mV/N]
hammer_sens_none = 2.361; % [mV/N]
sens_hammer = hammer_sens_steel; % default choice

sens_acc = 10.2; % [mV/(m/s^2)]
fmax_acc = 3000; % [Hz]

% definition of EU datas
for ii = 1:nset
    for jj = 1:ntest
        test{ii}.data{jj}.Dati(:,1) = test{ii}.data{jj}.Dati(:,1)./sens_hammer;
        test{ii}.data{jj}.Dati(:,2) = test{ii}.data{jj}.Dati(:,2)./sens_acc;
        test{ii}.data{jj}.Dati(:,3) = test{ii}.data{jj}.Dati(:,3)./sens_acc;
        test{ii}.data{jj}.Dati(:,4) = test{ii}.data{jj}.Dati(:,4)./sens_acc;
        test{ii}.data{jj}.Dati(:,5) = test{ii}.data{jj}.Dati(:,5)./sens_acc;
        test{ii}.data{jj}.Dati(:,6) = test{ii}.data{jj}.Dati(:,6)./sens_acc;
    end
end

% f sampling
safety_factor_aliasing = 2.56;
% fsamp = 810 * safety_factor_aliasing; % seen from mode 6 abaqus 
fsamp = 2500;

% creation of time vector
dt = 1/fsamp;
nsamples = length(test{1}.data{1}.Dati(:,1));
time = (0:dt:(nsamples*dt-dt))';

%% windowing of time domain signals - Tukey and Exponential
treshold = 0.02; % N
duration_test = 25; % seconds
pretrigger_duration = 0.1; % seconds

for ii = 1:nset
    for jj = 1:ntest
        idx = find(test{ii}.data{jj}.Dati(:,1) > treshold, 1); % index at which input signal is higher than treshold
        idx_list = (idx-pretrigger_duration*fsamp):1:(idx - pretrigger_duration*fsamp + duration_test*fsamp);
        test{ii}.data{jj}.Dati = [test{ii}.data{jj}.Dati(idx_list,1), test{ii}.data{jj}.Dati(idx_list,2), test{ii}.data{jj}.Dati(idx_list,3), test{ii}.data{jj}.Dati(idx_list,4), test{ii}.data{jj}.Dati(idx_list,5), test{ii}.data{jj}.Dati(idx_list,6)];
        % data{ii}.Dati = [data{ii}.Dati(idx_list,1), data{ii}.Dati(idx_list,2), data{ii}.Dati(idx_list,3)];
    end
end

N = length(test{1}.data{1}.Dati(: ,1));
% windowing with tukey window
w_tukey = [tukeywin(2*pretrigger_duration*fsamp); zeros((N-2*pretrigger_duration*fsamp), 1)];

% exponential window 
T = N/fsamp;          % seconds duration
t = linspace(0,(T-1), N);    % time axis
t0 = pretrigger_duration;         % initial time sample
P = 0.01;                       % percentage of final value
tau = (T - t0) / log(1/P);      % time constant

% defining exponential window 
w_exp = ones(size(t));
w_exp(t >= t0) = exp(-(t(t >= t0) - t0)/tau);
w_exp = w_exp';

% windowing 
for ii = 1:nset
    for jj = 1:ntest
        for kk = 1:n_acc
            if kk == 1
                test{ii}.data{jj}.Dati(:,kk) = w_tukey.*test{ii}.data{jj}.Dati(:,kk);
            else
                % data{ii}.Dati(:,jj) = w_exp.*data{ii}.Dati(:,jj);
            end
        end
    end
end

%% plotting time - domain

% set = 5; % default
% 
% t1 = linspace(0,(T), N);
% 
% plotting_title = [{'Input'}, {'Acc1'}, {'Acc2'}, {'Acc3'}, {'Acc4'}, {'Acc5'}];
% for jj = 1:ntest
%     figure('Name', ['Test ' num2str(jj)], 'NumberTitle', 'off');
%     plotting_cell = {test{set}.data{jj}.Dati(:, 1), test{set}.data{jj}.Dati(:, 2), test{set}.data{jj}.Dati(:, 3), test{set}.data{jj}.Dati(:, 4), test{set}.data{jj}.Dati(:, 5), test{set}.data{jj}.Dati(:, 6)};
%     for kk = 1:(n_acc + 1)
%         subplot(n_acc+1,1,kk)
%         plot(t1, plotting_cell{kk}, 'LineWidth', 1.2)
%         title(plotting_title{kk})
%         xlabel ('time [s]')
%         ylabel('EU')
%         xlim([0, 25])
%         grid on 
%     end
%     sgtitle(['Time domain plots: set ', num2str(set), 20]);
% end

%% Frequency domain calculation
% preallocation of set struct 
FFT = cell(nset, 1);
Pspectra = cell(nset, 1);
Cspectra = cell(nset, 1);
CC_Cspectra = cell(nset, 1);
Phase = cell(nset, 1);
H1 = cell(nset, 1);
H2 = cell(nset, 1);
Coherence = cell(nset, 1);
neg = 0; %to cancel out negative frequencies

for ii = 1:nset             % for each set 
    % preallocation of ffts window 
    FFT{ii}.input_ffts = cell(ntest,1);
    FFT{ii}.acc1_ffts = cell(ntest, 1);
    FFT{ii}.acc2_ffts = cell(ntest, 1);
    FFT{ii}.acc3_ffts = cell(ntest, 1);
    FFT{ii}.acc4_ffts = cell(ntest, 1);
    FFT{ii}.acc5_ffts = cell(ntest, 1);
    
    % population of the tests
    for jj = 1:ntest
        FFT{ii}.input_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 1));
        FFT{ii}.acc1_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 2));
        FFT{ii}.acc2_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 3));
        FFT{ii}.acc3_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 4));
        FFT{ii}.acc4_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 5));
        FFT{ii}.acc5_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 6));
    end
    
    % Power spectrum calculation
    Pspectra{ii}.input_Ps = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.input_ffts, neg);
    Pspectra{ii}.acc_Ps1 = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.acc1_ffts, neg);
    Pspectra{ii}.acc_Ps2 = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.acc2_ffts, neg);
    Pspectra{ii}.acc_Ps3 = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.acc3_ffts, neg);
    Pspectra{ii}.acc_Ps4 = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.acc4_ffts, neg);
    Pspectra{ii}.acc_Ps5 = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.acc5_ffts, neg);
    
    % Cross spectrum calculation
    Cspectra{ii}.acc_Cs1 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc1_ffts, neg);
    Cspectra{ii}.acc_Cs2 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc2_ffts, neg);
    Cspectra{ii}.acc_Cs3 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc3_ffts, neg);
    Cspectra{ii}.acc_Cs4 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc4_ffts, neg);
    Cspectra{ii}.acc_Cs5 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc5_ffts, neg);
    
    % complex conjugate cross spectrum 
    CC_Cspectra{ii}.acc_Cs1 = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs2 = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs3 = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs4 = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs5 = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.input_ffts, neg);
    
    % Phase 
    Phase{ii}.phase_acc1 = angle(Cspectra{ii}.acc_Cs1);
    Phase{ii}.phase_acc2 = angle(Cspectra{ii}.acc_Cs2);
    Phase{ii}.phase_acc3 = angle(Cspectra{ii}.acc_Cs3);
    Phase{ii}.phase_acc4 = angle(Cspectra{ii}.acc_Cs4);
    Phase{ii}.phase_acc5 = angle(Cspectra{ii}.acc_Cs5);
    
    
    % ------> Extimation of FRFs ------<
    
    % H1
    H1{ii}.abs_out1 = abs(Cspectra{ii}.acc_Cs1)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out2 = abs(Cspectra{ii}.acc_Cs2)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out3 = abs(Cspectra{ii}.acc_Cs3)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out4 = abs(Cspectra{ii}.acc_Cs4)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out5 = abs(Cspectra{ii}.acc_Cs5)./Pspectra{ii}.input_Ps;

    % H2
    H2{ii}.abs_out1 = abs(Pspectra{ii}.acc_Ps1)./abs(CC_Cspectra{ii}.acc_Cs1);
    H2{ii}.abs_out2 = abs(Pspectra{ii}.acc_Ps2)./abs(CC_Cspectra{ii}.acc_Cs2);
    H2{ii}.abs_out3 = abs(Pspectra{ii}.acc_Ps3)./abs(CC_Cspectra{ii}.acc_Cs3);
    H2{ii}.abs_out4 = abs(Pspectra{ii}.acc_Ps4)./abs(CC_Cspectra{ii}.acc_Cs4);
    H2{ii}.abs_out5 = abs(Pspectra{ii}.acc_Ps5)./abs(CC_Cspectra{ii}.acc_Cs5);
    
    % coherence 
    Coherence{ii}.gamma_acc1 = (Cspectra{ii}.acc_Cs1.*CC_Cspectra{ii}.acc_Cs1)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps1);
    Coherence{ii}.gamma_acc2 = (Cspectra{ii}.acc_Cs2.*CC_Cspectra{ii}.acc_Cs2)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps2);
    Coherence{ii}.gamma_acc3 = (Cspectra{ii}.acc_Cs3.*CC_Cspectra{ii}.acc_Cs3)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps3);
    Coherence{ii}.gamma_acc4 = (Cspectra{ii}.acc_Cs4.*CC_Cspectra{ii}.acc_Cs4)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps4);
    Coherence{ii}.gamma_acc5 = (Cspectra{ii}.acc_Cs5.*CC_Cspectra{ii}.acc_Cs5)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps5);
end

%% plotting H1, H2, Coherence and Phase for a chosen couple input-sensor
% frequency vector 
f = linspace(0, 900, length(Pspectra{1}.input_Ps)*(900/(fsamp/2))); % stop after six peaks
f2 = linspace(0, fsamp/2, length(Pspectra{1}.input_Ps)); % overall
% user input 
idx_set = input('Select a input (1-5): ');
idx_acc = input('Select a sensor (1-5): ');

% preallocating plotting vector
plotting_vector_h1 = {H1{idx_set}.abs_out1, H1{idx_set}.abs_out2, H1{idx_set}.abs_out3, H1{idx_set}.abs_out4, H1{idx_set}.abs_out5};
plotting_vector_h2 = {H2{idx_set}.abs_out1, H2{idx_set}.abs_out2, H2{idx_set}.abs_out3, H2{idx_set}.abs_out4, H2{idx_set}.abs_out5};
plotting_vector_coh = {Coherence{idx_set}.gamma_acc1, Coherence{idx_set}.gamma_acc2, Coherence{idx_set}.gamma_acc3, Coherence{idx_set}.gamma_acc4, Coherence{idx_set}.gamma_acc5};
plotting_vector_phase = {Phase{idx_set}.phase_acc1, Phase{idx_set}.phase_acc2, Phase{idx_set}.phase_acc3, Phase{idx_set}.phase_acc4, Phase{idx_set}.phase_acc5};

figure
subplot(2,1,1)
semilogy(f, plotting_vector_h1{idx_acc}(1:length(f)), 'LineWidth', 1.2,'Color', 'b')
hold on 
semilogy(f, plotting_vector_h2{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'r')
xlabel('frequency [Hz]')
ylabel('abs')

titl = title(['H1, H2 and Coherence for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)], 'FontSize', 13);
titl.Units = 'normalized';                         % Imposta unità in coordinate normalizzate
titl.Position(2) = titl.Position(2) + 0.05; 


yyaxis right
ax = gca;
ax.YAxis(2).Color = 'k';
plot(f,plotting_vector_coh{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'k')
ylim([0, 3])
yline(1, '--k', 'LineWidth', 1.2)
grid on 
legend('H1', 'H2', 'Coherence', 'Location','northwest')

ylabel('Coherence')
hold off

subplot(2,1,2)
plot(f,plotting_vector_phase{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'b')
grid on
xlabel('frequency [Hz]')
ylabel('angle [rad]')
% title(['Phase for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)], 'FontSize', 15)

%%
% 
% % frequency vector 
% f = linspace(0, 900, length(Pspectra{1}.input_Ps)*(900/(fsamp/2))); % stop after six peaks
% f2 = linspace(0, fsamp/2, length(Pspectra{1}.input_Ps)); % overall
% % user input 
% idx_set = input('Select a input (1-5): ');
% idx_acc = input('Select a sensor (1-5): ');
% 
% % preallocating plotting vector
% plotting_vector_h1 = {H1{idx_set}.abs_out1, H1{idx_set}.abs_out2, H1{idx_set}.abs_out3, H1{idx_set}.abs_out4, H1{idx_set}.abs_out5};
% plotting_vector_h2 = {H2{idx_set}.abs_out1, H2{idx_set}.abs_out2, H2{idx_set}.abs_out3, H2{idx_set}.abs_out4, H2{idx_set}.abs_out5};
% plotting_vector_coh = {Coherence{idx_set}.gamma_acc1, Coherence{idx_set}.gamma_acc2, Coherence{idx_set}.gamma_acc3, Coherence{idx_set}.gamma_acc4, Coherence{idx_set}.gamma_acc5};
% plotting_vector_phase = {Phase{idx_set}.phase_acc1, Phase{idx_set}.phase_acc2, Phase{idx_set}.phase_acc3, Phase{idx_set}.phase_acc4, Phase{idx_set}.phase_acc5};
% 
% figure
% semilogy(f, plotting_vector_h1{idx_acc}(1:length(f)), 'LineWidth', 1.2,'Color', 'b')
% hold on 
% semilogy(f, plotting_vector_h2{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'r')
% xlabel('frequency [Hz]')
% ylabel('abs')
% 
% titl = title(['H1, H2 and Coherence for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)], 'FontSize', 13);
% titl.Units = 'normalized';                         % Imposta unità in coordinate normalizzate
% titl.Position(2) = titl.Position(2) + 0.05; 
% legend('H1', 'H2', 'Location','northwest')
% 
% %xlim([120, 220])




%% plotting power and cross spectrum
% % power spectra plot -- OPTIONAL 
% plotting_cell = {input_Ps, acc_Ps1, acc_Ps2, acc_Ps3, acc_Ps4, acc_Ps5};
% % plotting_cell = {input_Ps, acc_Ps1, acc_Ps2};
% plotting_title = [{'PS-input'}, {'PS-out1'}, {'PS-out2'}, {'PS-out3'}, {'PS-out4'}, {'PS-out5'}];
% 
% figure('Name', 'Power Spectra', 'NumberTitle', 'off');
% for jj = 1:(n_acc+1)
%     subplot(n_acc+1,1,jj)
%     semilogy(f,abs(plotting_cell{jj}), 'LineWidth', 1.2)
%     grid on
%     hold on 
%     xlabel('f')
%     ylabel('EU')
%     title(plotting_title(jj))
% end
% 
% % cross spectra plot
% % plotting_cell = {input_Ps, acc_Ps1, acc_Ps2, acc_Ps3, acc_Ps4, acc_Ps5};
% plotting_cell = {input_Ps, acc_Cs1, acc_Cs2, acc_Cs3, acc_Cs4, acc_Cs5};
% plotting_title = [{'PS-input'}, {'CS-out1'}, {'CS-out2'}, {'CS-out3'}, {'CS-out4'}, {'CS-out5'}];
% 
% figure('Name', 'Cross Spectra', 'NumberTitle', 'off');
% for jj = 1:(n_acc+1)
%     subplot(n_acc+1,1,jj)
%     semilogy(f,abs(plotting_cell{jj}), 'LineWidth', 1.2)
%     grid on
%     hold on 
%     xlabel('f')
%     ylabel('EU')
%     title(plotting_title(jj))
% end


%% MDOF 

% We choose H1 because we believe we have more noise on the output
n_inp = n_acc;
Pspectra_n = cell(nset, 1);
Cspectra_n = cell(nset, 1);
neg = 1; %now we want to keep negative frequencies for doing the inverse F.T.
H1_C = cell(n_inp, 1);

for ii = 1:nset 
    % Power spectrum calculation
    Pspectra_n{ii}.input_Ps_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.input_ffts, neg);
    Pspectra_n{ii}.acc_Ps1_n = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.acc1_ffts, neg);
    Pspectra_n{ii}.acc_Ps2_n = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.acc2_ffts, neg);
    Pspectra_n{ii}.acc_Ps3_n = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.acc3_ffts, neg);
    Pspectra_n{ii}.acc_Ps4_n = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.acc4_ffts, neg);
    Pspectra_n{ii}.acc_Ps5_n = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.acc5_ffts, neg);

    % Cross spectrum calculation
    Cspectra_n{ii}.acc_Cs1_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc1_ffts, neg);
    Cspectra_n{ii}.acc_Cs2_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc2_ffts, neg);
    Cspectra_n{ii}.acc_Cs3_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc3_ffts, neg);
    Cspectra_n{ii}.acc_Cs4_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc4_ffts, neg);
    Cspectra_n{ii}.acc_Cs5_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc5_ffts, neg);

    % H1-COMPLEX
    H1_C{ii, 1}.out1 = Cspectra_n{ii}.acc_Cs1_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out2 = Cspectra_n{ii}.acc_Cs2_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out3 = Cspectra_n{ii}.acc_Cs3_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out4 = Cspectra_n{ii}.acc_Cs4_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out5 = Cspectra_n{ii}.acc_Cs5_n./Pspectra_n{ii}.input_Ps_n;

end

%% Matrix h for each time instant

h = zeros(n_inp, n_acc, N);

for ii = 1:n_inp
    h(1, ii, :) = ifft(H1_C{ii, 1}.out1);
    h(2, ii, :) = ifft(H1_C{ii, 1}.out2);
    h(3, ii, :) = ifft(H1_C{ii, 1}.out3);
    h(4, ii, :) = ifft(H1_C{ii, 1}.out4);
    h(5, ii, :) = ifft(H1_C{ii, 1}.out5);
end

%%

% % if you want to see the unit impulse response 
% figure
% B = squeeze(h(1,5,:));
% plot(t, B)

%% Now we cut the h at t = 20sec

t_cut = 0:dt:((20*length(t)/(25))*dt-dt);
h_cut = h(:, :, 1:(20*length(t)/(25)));
figure
B_cut = squeeze(h_cut(2,5,:));
plot(t_cut, B_cut)

%% Now we define Hmn

% NB!!!!!!!!!!!
%THIS PART HAS NOT BEEN CANCELLED BECAUSE I'M NOT SURE, SO I KEEP IT BUT
%THE PART THAT WE HAVE TO RUN IS Whole stabilization diagram


%tot = length(t_cut);
% m = 1000;
% n = tot - m;
% 
% Hmn = zeros(m*n_acc, n*n_inp);
% 
% for ii = 1:m
%     for jj = 1:n
%         Hmn((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2);
%     end
% end
% 
% 
% 
% 
% Hmn_next= zeros(m*n_acc, n*n_inp);
% 
% for ii = 1:m
%     for jj = 1:n
%         Hmn_next((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2 + 1);
%     end
% end

%% Now we do the diagonalization


% NB!!!!!!!!!!!
%THIS PART HAS NOT BEEN CANCELLED BECAUSE I'M NOT SURE, SO I KEEP IT BUT
%THE PART THAT WE HAVE TO RUN IS Whole stabilization diagram


% Hmn_pseudo = pinv(Hmn);
% Hmm_final = Hmn_next*Hmn_pseudo;
% 
% [phi, eigenvalues] = eig(Hmm_final);
% 
% poles = (log(diag(eigenvalues)))./dt;
% freq = abs(poles)./(2*pi);
% 
% 
% x = 1;                 % larghezza desiderata per ogni rettangolo
% 
% % Calcolo i bin edges con passo x
% edges = min(freq):x:810;  % +x per includere l'ultimo dato
% 
% figure
% % Crea l'istogramma
% histogram(freq, 'BinEdges', edges);


%% Frequency stabilization:

tot = length(t_cut);
m_values = [12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132];
% m_values = [12, 84, 156, 228, 300, 400];
freq_stab = [];
poles_tot = [];
phi = [];

for idx = 1:length(m_values)
    m = m_values(idx);
    n = tot - m;

    Hmn = zeros(m*n_acc, n*n_inp);
    Hmn_next = zeros(m*n_acc, n*n_inp);

    for ii = 1:m
        for jj = 1:n
            Hmn((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2);
            Hmn_next((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2 + 1);
        end
    end

    H_pinv = pinv(Hmn);
    Hmm = Hmn_next * H_pinv;
    [eigenvectors, eigenvalues] = eig(Hmm);
    poles = (log(diag(eigenvalues)))./dt;
    freq = abs(poles)./(2*pi);
    poles_tot = [poles_tot; poles];
    freq_stab = [freq_stab; freq];
    phi_i = eigenvectors(1:5, :);
    phi = [phi, phi_i];
end

freq_res = 0.6; %width of the frequency range for the stabilization diagram 0.1%
edges = min(freq):freq_res:820;

figure
histogram(freq_stab, 'BinEdges', edges)


%% MDOF - (Previous section, kept for context)
% ... (previous code for H1_C, h, h_cut) ...

% Matrix h for each time instant (already correctly calculated)
% h = zeros(n_inp, n_acc, N);
% for ii = 1:n_inp
%     h(1, ii, :) = ifft(H1_C{ii, 1}.out1);
%     h(2, ii, :) = ifft(H1_C{ii, 1}.out2);
%     h(3, ii, :) = ifft(H1_C{ii, 1}.out3);
%     h(4, ii, :) = ifft(H1_C{ii, 1}.out4);
%     h(5, ii, :) = ifft(H1_C{ii, 1}.out5);
% end

% Now we cut the h at t = 20sec (already correctly calculated)
% t_cut = 0:dt:((20*length(t)/(25))*dt-dt);
% h_cut = h(:, :, 1:(20*length(t)/(25)));

% --- Assicurati che queste variabili siano definite prima di questo blocco ---
% h_cut: matrice delle risposte all'impulso tagliate (n_inp x n_acc x N_samples)
% dt: passo temporale
% n_acc: numero di accelerometri (output)
% n_inp: numero di input (assumendo 1 per il martello, come discusso)
% -------------------------------------------------------------------------

%% Modal Identification using ITD (Ibrahim Time Domain) with SVD for Reduction

% % Adattamento di n_inp a 1 se l'input è un singolo martello
% % Se 'h_cut' è stata generata con un solo input, la sua prima dimensione sarà 1.
% if size(h_cut, 1) ~= n_inp
%     warning('Mismatch between actual h_cut dimensions and n_inp. Setting n_inp to size(h_cut, 1).');
%     n_inp = size(h_cut, 1); 
% end
% 
% % Flattening delle risposte all'impulso per la costruzione della matrice di Hankel.
% % h_flat avrà dimensioni: (n_inp * n_acc) x N_samples_h_cut
% h_flat = reshape(h_cut, n_inp * n_acc, size(h_cut, 3)); 
% 
% N_samples_h_cut = size(h_cut, 3); % Numero totale di campioni temporali in h_cut
% 
% % Definizione della profondità del blocco della matrice di Hankel (r).
% % Questo 'r' dovrebbe essere sufficientemente grande per catturare la dinamica
% % del sistema, spesso da 2 a 5 volte il massimo ordine del modello atteso.
% % Per iniziare, usiamo un valore che tenga conto di un ampio range di ordini del modello.
% r = 250; % Puoi aumentare questo valore (es. 500, 1000) se i tuoi dati lo permettono
%          % e se vedi che modi di ordine superiore non sono catturati.
% 
% % Controllo per assicurarsi che 'r' non superi la lunghezza dei dati
% if r >= N_samples_h_cut
%     error('Il valore di "r" (profondità blocco Hankel) è troppo grande. Riduci "r" o aumenta la durata del segnale in h_cut.');
% end
% 
% num_data_channels = n_inp * n_acc; % Numero totale di "canali" di risposta (es. 1x5=5 per singolo input/multi-output)
% num_columns_hankel_matrix = N_samples_h_cut - r; % Numero di colonne nella matrice di osservazione
% 
% % Controllo per assicurarsi che ci siano abbastanza colonne
% if num_columns_hankel_matrix <= 0
%     error('Numero insufficiente di campioni in h_cut per costruire la matrice di Hankel con il "r" specificato.');
% end
% 
% % Preallocazione delle matrici di osservazione O_0 e O_1 (shifted)
% % O_0: (num_data_channels * r) x num_columns_hankel_matrix
% % O_1: (num_data_channels * r) x num_columns_hankel_matrix
% O_0 = zeros(num_data_channels * r, num_columns_hankel_matrix);
% O_1 = zeros(num_data_channels * r, num_columns_hankel_matrix);
% 
% % Popolamento delle matrici di osservazione
% fprintf('Costruzione delle matrici di Hankel O_0 e O_1...\n');
% for j_col = 1:num_columns_hankel_matrix
%     % Blocco corrente di risposte per O_0
%     current_block = h_flat(:, j_col : (j_col + r - 1));
%     O_0(:, j_col) = current_block(:); % Appiattisci in una colonna
% 
%     % Blocco shiftato di risposte per O_1
%     shifted_block = h_flat(:, (j_col + 1) : (j_col + r));
%     O_1(:, j_col) = shifted_block(:); % Appiattisci in una colonna
% end
% 
% % Applicazione della SVD alla matrice di osservazione O_0
% fprintf('Esecuzione della SVD su O_0 (matrice di osservazione)...\n');
% [U_svd, S_svd, ~] = svd(O_0, 'econ'); % 'econ' per SVD in formato economico
% 
% % Definizione degli ordini del modello da testare per lo stabilization diagram
% % Inizia da un valore basso e arriva a un valore alto.
% % Scegli un range adeguato in base alla complessità attesa del sistema
% % e al numero massimo di modi che pensi di poter identificare.
% model_orders = 2:2:min(size(U_svd, 2), 2 * 300); % Ordini pari fino a 600 (o al rango disponibile)
% % Puoi aumentare il 300 a un valore maggiore se hai molti più modi o un sistema complesso.
% 
% % Preallocazione degli array per salvare tutti i poli identificati
% all_identified_poles_freq = [];
% all_identified_poles_damping = [];
% all_identified_poles_order = [];
% 
% fprintf('Calcolo dei poli per diversi ordini del modello (per lo Stabilization Diagram)...\n');
% for p_idx = 1:length(model_orders)
%     p = model_orders(p_idx); % Ordine del modello corrente
% 
%     % Assicurati che 'p' non superi il rango effettivo dalla SVD
%     if p > size(U_svd, 2)
%         warning('Ordine del modello %d supera il numero di valori singolari disponibili. Interruzione.', p);
%         break; % Esci dal ciclo se superi il rango massimo identificabile
%     end
% 
%     % Estrazione dei blocchi di U per il calcolo della matrice di stato A_sys
%     % U_p_block1 corrisponde alla parte superiore della matrice di osservabilità troncata
%     % U_p_block2 corrisponde alla parte inferiore/shiftata della matrice di osservabilità troncata
% 
%     % Check per dimensioni valide prima di estrarre i blocchi
%     if (num_data_channels * (r - 1)) <= 0 || (num_data_channels + 1) > size(U_svd, 1)
%         warning('Dimensioni non valide per l''estrazione dei blocchi U_p_block per l''ordine %d. Controllare r o num_data_channels. Salto questo ordine.', p);
%         continue; 
%     end
% 
%     U_p_block1 = U_svd(1:(num_data_channels * (r - 1)), 1:p); 
%     U_p_block2 = U_svd((num_data_channels + 1):end, 1:p); 
% 
%     % Verifica che le dimensioni delle colonne siano coerenti
%     if size(U_p_block1, 2) ~= size(U_p_block2, 2)
%         warning('Discrepanza di dimensione nelle colonne dei blocchi U per l''ordine %d. Salto questo ordine.', p);
%         continue;
%     end
% 
%     % Calcolo della matrice di transizione di stato A_sys
%     % A_sys ha dimensione p x p
%     A_sys = U_p_block2 * pinv(U_p_block1);
% 
%     % Verifica la validità di A_sys (può essere NaN o Inf per problemi numerici)
%     if isempty(A_sys) || any(isnan(A_sys(:))) || any(isinf(A_sys(:)))
%         warning('Matrice di stato A_sys non valida per l''ordine %d. Salto questo ordine.', p);
%         continue;
%     end
% 
%     % Calcolo autovalori (e quindi i poli) di A_sys
%     [~, eigenvalues] = eig(A_sys);
% 
%     % Estrazione dei poli (frequenze complesse)
%     poles_current_order = (log(diag(eigenvalues)))./dt;
% 
%     % Filtra per mantenere solo un polo per ogni coppia complessa coniugata
%     % (solitamente si prendono quelli con parte immaginaria non negativa)
%     valid_poles_idx_current = imag(poles_current_order) >= 0;
%     poles_current_order = poles_current_order(valid_poles_idx_current);
% 
%     % Calcolo frequenza e smorzamento dai poli
%     freq_current = abs(poles_current_order) ./ (2*pi);
%     damping_ratio_current = -real(poles_current_order) ./ abs(poles_current_order);
% 
%     % Aggiungi tutti i poli validi (anche quelli "numerici" ma non NaN/Inf)
%     % per l'ordine corrente agli array generali.
%     % Escludiamo solo i poli che sono NaN o Inf, poiché causerebbero problemi di plotting.
%     valid_plot_idx = ~isnan(freq_current) & ~isinf(freq_current) & ...
%                      ~isnan(damping_ratio_current) & ~isinf(damping_ratio_current);
% 
%     all_identified_poles_freq = [all_identified_poles_freq; freq_current(valid_plot_idx)];
%     all_identified_poles_damping = [all_identified_poles_damping; damping_ratio_current(valid_plot_idx)];
%     all_identified_poles_order = [all_identified_poles_order; repmat(p, sum(valid_plot_idx), 1)];
% 
%     fprintf('Identificati %d poli (filtrati per validità plot) per l''ordine %d.\n', sum(valid_plot_idx), p);
% end
% 
% fprintf('Calcolo dei poli completato. Inizio creazione diagramma.\n');
% 
% % --- Creazione del Plot del Diagramma di Stabilizzazione ---
% % Questo plot mostra tutti i poli identificati, senza selezione di stabilità.
% % Saranno la base per la tua selezione manuale o automatica.
% 
% figure('Name', 'Diagramma di Stabilizzazione (Frequenza vs. Ordine)', 'NumberTitle', 'off');
% scatter(all_identified_poles_order, all_identified_poles_freq, 20, 'b', 'filled', 'o', 'MarkerFaceAlpha', 0.5);
% xlabel('Ordine del Modello');
% ylabel('Frequenza [Hz]');
% title('Diagramma di Stabilizzazione (Tutti i Poli Identificati)');
% grid on;
% % Puoi impostare limiti Y se conosci un range di frequenze di interesse
% % Es: ylim([0 1000]); 
% set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
% 
% figure('Name', 'Diagramma di Stabilizzazione (Frequenza vs. Smorzamento)', 'NumberTitle', 'off');
% scatter(all_identified_poles_freq, all_identified_poles_damping, 20, 'r', 'filled', 'o', 'MarkerFaceAlpha', 0.5);
% xlabel('Frequenza [Hz]');
% ylabel('Rapporto di Smorzamento');
% title('Diagramma di Stabilizzazione (Frequenza vs. Smorzamento)');
% grid on;
% % Puoi impostare limiti X e Y se conosci un range di frequenze/smorzamenti di interesse
% % Es: xlim([0 1000]); ylim([0 0.1]); 
% set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
% 
% fprintf('I dati dei poli sono disponibili in:\n');
% fprintf('  - all_identified_poles_freq\n');
% fprintf('  - all_identified_poles_damping\n');
% fprintf('  - all_identified_poles_order\n');
% fprintf('Pronti per la tua successiva fase di selezione dei poli stabili.\n');
%% Damping stabilization:
freq_stable_poles = cell(6, 1); %we put here all the poles that are stable in frequency
freq_stable_modes = cell(6, 1);
damping_cell = cell(6, 1);
freq_range = [247.2, 327.6, 391.8, 596.4, 620.4, 806.4];


for ii = 1:length(freq_range)
    kk = 0;
    for jj = 1:length(poles_tot)
        
        current_mode = phi(:, jj);
        current_pole = poles_tot(jj);
        current_freq = abs(current_pole)./(2*pi);

        if(current_freq >= freq_range(ii) && current_freq <= (freq_range(ii)+freq_res))
            kk = kk + 1;
            freq_stable_poles{ii, 1}(kk) = current_pole;
            freq_stable_modes{ii, 1}(:, kk) = current_mode;
            %we compile the damping cell 
            current_damping = (-1)*real(current_pole)./(current_freq*2*pi);
            damping_cell{ii, 1}(kk) = current_damping;
        end
    end
end


mostFrequentInterval_damping = zeros(6, 2);
%now we plot the stabilization diagrams for each mode for the damping
for ii = 1:6
    max_damping = max(damping_cell{ii, 1}(:));
    damping_res = (0.5/100)*max_damping;
    edges = min(damping_cell{ii, 1}(:)) : damping_res : max_damping;

    figure
    histogram(damping_cell{ii, 1}(:), 'BinEdges', edges)

    % Calcolo dell'istogramma (senza plottarlo)
    [counts, ~] = histcounts(damping_cell{ii, 1}(:), edges);

    % Trovo il massimo conteggio
    [maxCount, idxMax] = max(counts);

    % Estrazione dell'intervallo associato al valore più frequente
    mostFrequentInterval_damping(ii,:) = [edges(idxMax), edges(idxMax+1)];
    
end

freq_damping_stable_poles = cell(6, 1);
freq_damping_stable_modes = cell(6, 1);

for ii = 1:6
    kk = 0;
    for jj = 1:length(freq_stable_poles{ii, 1}(:))

        current_mode = freq_stable_modes{ii, 1}(:, jj);
        current_freq = abs(freq_stable_poles{ii, 1}(jj))./(2*pi);
        current_damping = (-1)*real(freq_stable_poles{ii, 1}(jj))./(current_freq*2*pi);

        if(current_damping >= mostFrequentInterval_damping(ii,1) && current_damping <= mostFrequentInterval_damping(ii,2))
            kk = kk + 1;

            freq_damping_stable_modes{ii, 1}(:, kk) = freq_stable_modes{ii, 1}(:, jj);
            freq_damping_stable_poles{ii, 1}(kk) = freq_stable_poles{ii, 1}(jj);
        end
    end
end


%% Modes

%siamo in caso di damping proporzionale perchè abbiamo un oggetto metallico
%e quindi ci aspettiamo di avere una parte immaginaria bassa nei modi,
%infatti è così!! Quindi quello che facciamo è scartare la parte
%immaginaria (un ordine di grandezza in meno) e tenere solo quella reale.
% Se avessimo avuto modi complessi con parte immaginaria non trascurabile 
%avremmo trovato una massa modale complessa e il rapporto tornava sensato
%di nuovo 

freq_damping_stable_modes_real = cell(6, 1);
for ii = 1:6
    freq_damping_stable_modes_real{ii, 1}(:, :) = real(freq_damping_stable_modes{ii, 1}(:, :));
end


MAC = [];
for ii = 1:6
    %we choose as a reference one mode and we compute the mac with respect
    %to this reference mode for all the other modes of the ii-th component
    reference_mode = freq_damping_stable_modes_real{ii, 1}(:, 1);
    k = 0;
    
    numerator = [];
    denominator = [];
    for jj = 1:length(freq_damping_stable_modes_real{ii, 1}(1, :))
        current_mode = freq_damping_stable_modes_real{ii, 1}(:, jj);

        numerator = abs(reference_mode' * current_mode)^2;
        denominator = (reference_mode' * reference_mode) * (current_mode' * current_mode);
        MAC(ii, jj) = numerator / denominator;


        if(MAC > 0.9)
            k = k + 1;
        end
    end
end
        

%SINCE WE DID A GREAT JOB BEFORE THE freq_damping_stable_poles ARE EXACTLY
%EQUAL TO THE POLES ASSOCIATED WITH THE POLES THAT HAVE STABLE MODES (MAC
%BIG ENOUGH FOR ALL) IL MAC NON SI ABBASSA MAI, RIMANE SEMPRE ALTO PERCHE
%ABBIAMO FILTRATO MOLTO NEI PASSAGGI PRECEDENTI DI FREQUENZA E DAMPING.

final_stable_poles = freq_damping_stable_poles;
final_stable_modes = freq_damping_stable_modes_real;


%% Averages:

f0_average = [];
csi_average = [];
modes_average = [];
poles_average = [];


for ii = 1:6
    
    %we have complex and conjugated poles, in order to do the average we
    %need to split in 2 sets of poles otherwise we get a real number
    %because the imaginary parts cancel out each other
    current_pole1 = final_stable_poles{ii, 1}(1:2:end);
    current_pole2 = final_stable_poles{ii, 1}(2:2:end);
    poles_average(2*ii - 1,1) = mean(current_pole1);
    poles_average(2*ii,1) = mean(current_pole2);


    current_f0 = abs(final_stable_poles{ii, 1}(:))./(2*pi);
    f0_average(2*ii - 1, 1) = mean(current_f0);
    f0_average(2*ii, 1) = mean(current_f0);

    current_csi = (-1)*real(final_stable_poles{ii, 1}(:))./(current_f0*2*pi);
    csi_average(2*ii-1, 1) = mean(current_csi);
    csi_average(2*ii, 1) = mean(current_csi);

    current_mode = final_stable_modes{ii, 1}(:, :);
    modes_average(2*ii - 1, :) = mean(current_mode, 2);
    modes_average(2*ii, :) = mean(current_mode, 2);

end



%% Modal mass:

N_modes = 6;
f_new = f(4:end);
A = zeros(length(f_new),(2*N_modes+2));
X = zeros((2*N_modes + 2), 1);
B = zeros(N_modes,1);

B = H1_C{3, 1}.out3(1:length(f_new)); %this is Hpq where p = sensor and q = input

for ii = 1:(2*N_modes)

    %in this case we use 3 because the sensor is the 3-th
    column_numerator = ones(length(f_new), 1).*modes_average(ii, 3).^2;
    column_denominator = 1i.*(2*pi.*f_new') - (ones(length(f_new), 1).*poles_average(ii));

    A(:, ii) = column_numerator./column_denominator;

end

A(:, end-1) = ones(length(f_new), 1);

A(:, end) = 1./((f_new.*2*pi).^2);


X = pinv(A) * B;  % pseudoinversa più robusta

X_2 = A\B;

%% Finding m_modal:

Q_vector = X(1:2:12, 1);
Q_vector_2 = X_2(1:2:12, 1);
modal_masses = [];
modal_masses_2 = [];

for ii = 1 : N_modes

    modal_masses(ii) = (-1*(2*pi.*f0_average(ii*2-1)).^2)/(Q_vector(ii)*2*1i*f0_average(ii*2-1)*2*pi*sqrt(1-(csi_average(ii*2-1))^2));
    modal_masses_2(ii) = (-1*(2*pi.*f0_average(ii*2-1)).^2)/(Q_vector_2(ii)*2*1i*f0_average(ii*2-1)*2*pi*sqrt(1-(csi_average(ii*2-1))^2));

end


%% FRF estimation:

H_pred_modal = 0;
H_pred_modal_2 = 0;

for ii = 1:6
    H_pred_i = (modes_average(2*ii, 3)*modes_average(2*ii, 3)*(2*pi.*f).^2)./(abs(modal_masses(ii))*(-(2*pi.*f).^2 + 2*1i.*(2*pi.*f).*(2*pi*f0_average(2*ii)).*csi_average(2*ii) + (2*pi*f0_average(2*ii))^2));
    H_pred_modal = H_pred_modal + H_pred_i;
    H_pred_i_2 = (modes_average(2*ii, 3)*modes_average(2*ii, 3)*((2*pi.*f).^2))./(abs(modal_masses_2(ii))*(-(2*pi.*f).^2 + 2*1i.*(2*pi.*f).*(2*pi*f0_average(2*ii)).*csi_average(2*ii) + (2*pi*f0_average(2*ii))^2));
    H_pred_modal_2 = H_pred_modal_2 + H_pred_i;
end

H_pred = A*X;

% Funzione di trasferimento sperimentale
H_exp = H1_C{3, 1}.out3(1:length(f_new));

% Plot combinato modulo + fase
figure

% --- Modulo (scala logaritmica) ---
subplot(2,1,1)
semilogy(f_new, abs(H_pred_modal(1:length(f_new))), 'r', 'LineWidth', 1.2);
hold on
semilogy(f_new, abs(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2);
semilogy(f_new, abs(H_pred(1:length(f_new))), 'LineWidth', 1.2);
semilogy(f_new, abs(H_exp), 'b', 'LineWidth', 1.2);
ylabel('|H(f)|')
title('Modulo e fase della funzione di trasferimento')
legend('Mdof pinv1', 'Mdof pinv2', 'Mdof pred', 'Sperimentale', 'Location', 'northeast')
grid on

% --- Fase (in gradi) ---
subplot(2,1,2)
plot(f_new, angle(H_pred_modal(1:length(f_new))), 'r', 'LineWidth', 1.2);
hold on
%plot(f_new, angle(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2);
%plot(f_new, angle(H_pred(1:length(f_new))), 'LineWidth', 1.2);
plot(f_new, angle(H_exp), 'b', 'LineWidth', 1.2);
xlabel('Frequenza [Hz]')
ylabel('Fase [°]')
legend('Mdof pinv1', 'Mdof pinv2', 'Mdof pred', 'Sperimentale', 'Location', 'northeast')
grid on

%% MAC final computation

idx = [12110, 11075, 7884, 10814, 11299];

fileNames = [{'Mode1.txt'},{'Mode2.txt'},{'Mode3.txt'},{'Mode4.txt'},{'Mode5.txt'},{'Mode6.txt'}];

N = 6; %number of investigated modes

InputMatrix = [];
SDOF_struct = [];
MDOF_modes = [];
SDOF_modes = [];
Norm_Experimental_modes = [];
Norm_MDOF_modes = [];
Norm_SDOF_modes = [];


for ii=1:N
    A = readmatrix(fileNames{ii});   %create matrix from txt file
    B = A(:,4);                      % extract last column (z axis)
    InputMatrix = [InputMatrix , B];   % populate input matrix
end

Experimental_modes = InputMatrix(idx,:);
MDOF_modes = modes_average(1:2:end, :)';
SDOF_struct = load('SDOF_modes.mat');
SDOF_modes = SDOF_struct.exp_modes_shape;

for ii = 1 : N_modes
    Norm_Experimental_modes(:, ii) = Experimental_modes(:, ii)./max(abs(Experimental_modes(:, ii)));
    Norm_MDOF_modes(:, ii) = MDOF_modes(:, ii)./max(abs(MDOF_modes(:, ii)));
    Norm_SDOF_modes(:, ii) = SDOF_modes(:, ii)./max(abs(SDOF_modes(:, ii)));

end

final_MAC_MDOF = [];
final_MAC_SDOF = [];

for ii = 1:N_modes
    numerator_MDOF = [];
    numerator_SDOF = [];
    denominator_MDOF = [];
    denominator_SDOF = [];

    numerator_MDOF = abs(Experimental_modes(:, ii)' * MDOF_modes(:, ii))^2;
    numerator_SDOF = abs(Experimental_modes(:, ii)' * SDOF_modes(:, ii))^2;
    denominator_MDOF = (Experimental_modes(:, ii)' * Experimental_modes(:, ii)) * (MDOF_modes(:, ii)' * MDOF_modes(:, ii));
    denominator_SDOF = (Experimental_modes(:, ii)' * Experimental_modes(:, ii)) * (SDOF_modes(:, ii)' * SDOF_modes(:, ii));

    final_MAC_MDOF(ii, :) = numerator_MDOF / denominator_MDOF;
    final_MAC_SDOF(ii, :) = numerator_SDOF / denominator_SDOF;
end



%%

load('FRF_Final_SDOF.mat')
load('FRF_Phase_SDOF.mat')
load('FRF_Final_MDOF_Prony.mat')

% H_pred_column = H_pred_modal';
% leng_SDOF = size(frf_formula);
% leng_MDOF = size(H_pred_column);
% 
% % slicing of vectors
% [~, idx_min_frf] = min(abs(f - 240));
% [~, idx_max_frf] = min(abs(f - 840));
% 
% frf_cutted = frf_approx(idx_min_frf:idx_max_frf);
% omega_cutted = omega(idx_min_frf:idx_max_frf);

% 
% figure;
% subplot(2,1,1)
% semilogy(abs(H_pred_column((leng_MDOF-leng_SDOF-1500):end)))
% hold on;
% semilogy(abs(frf_formula))
% semilogy(abs(H_exp((leng_MDOF-leng_SDOF-1500):end)))
% legend('MDOF', 'SDOF', 'Location', 'northeast')
% 
% subplot(2,1,2)
% plot(-angle(H_pred_column((leng_MDOF-leng_SDOF-1500):end)))
% hold on;
% plot(angle(frf_approx((leng_MDOF-leng_SDOF-1500):end))+pi)
% plot(unwrap(angle(H_exp((leng_MDOF-leng_SDOF-1500):end))))
% 
% legend('MDOF', 'SDOF', 'Location', 'northeast')

% --- Assicurati che queste variabili siano già caricate e disponibili ---
% f: Vettore delle frequenze (globale, prima del taglio)
% omega: Vettore delle frequenze angolari (globale, prima del taglio)
% H_pred_modal: La tua FRF predetta dal modello MDOF
% frf_formula: La tua FRF dal modello SDOF (Modulo)
% frf_approx: La tua FRF dal modello SDOF (Fase)
% H1_C: Cella contenente la FRF sperimentale (es. H1_C{3,1}.out3)
% ----------------------------------------------------------------------

% Caricamento dati (se non già caricati prima)
% load('FRF_Final_SDOF.mat')   % Assicurati che contenga frf_formula e f, omega
% load('FRF_Phase_SDOF.mat')   % Assicurati che contenga frf_approx
% Assumo che H_pred_modal e H1_C siano caricate o generate prima

% Definizione del range di frequenze per il taglio
% Questi valori (240 e 840 Hz) sono i tuoi punti di interesse
start_freq_cut = 240;
end_freq_cut = 840;

% Calcolo degli indici per il taglio basato sul vettore delle frequenze 'f'
[~, idx_min_frf] = min(abs(f - start_freq_cut));
[~, idx_max_frf] = min(abs(f - end_freq_cut));

omega = f.*2*pi;
frf_cutted = frf_approx(idx_min_frf:idx_max_frf);
omega_cutted = omega(idx_min_frf:idx_max_frf);
frf_prony_cutted = H_pred_modal_prony';

% Adattamento H_pred_modal e estrazione H_exp
H_pred_column = H_pred_modal'; % Assicurati che sia un vettore colonna se necessario

% Esecuzione del taglio su tutti i vettori rilevanti
frf_approx_cutted = frf_approx(idx_min_frf:idx_max_frf);
H_pred_column_cutted = H_pred_column(idx_min_frf:idx_max_frf);
H_exp_cutted = H1_C{3, 1}.out3(idx_min_frf:idx_max_frf); % Preleva la FRF sperimentale tagliata
f_cutted = f(idx_min_frf:idx_max_frf); % Anche il vettore frequenze va tagliato!

H_exp_cutted_phase = plotting_vector_phase{1, 3}(idx_min_frf:idx_max_frf);

% % --- INIZIO CODICE PLOT MIGLIORATO ---
% 
% figure('Name', 'Confronto FRF Modulo e Fase - Modelli', 'NumberTitle', 'off'); % Titolo generale della finestra
% 
% % --- Subplot 1: Modulo (scala logaritmica) ---
% subplot(2,1,1);
% semilogy(f_cutted, abs(H_pred_column_cutted), 'DisplayName', 'MDOF', 'LineWidth', 1.5, 'Color', 'r'); % Blu MATLAB
% hold on;
% % semilogy(f_cutted, abs(frf_formula), 'DisplayName', 'SDOF', 'LineWidth', 1.5, 'Color', 'r'); % Arancio MATLAB
% semilogy(f_cutted, abs(H_exp_cutted), 'DisplayName', 'Experimental', 'LineWidth', 1.5, 'Color', 'b'); % Verde MATLAB
% 
% xlabel('Frequenza [Hz]', 'FontSize', 11);
% ylabel('Modulo |H(f)|', 'FontSize', 11);
% title('Module FRF Comparison: Models vs. Experimental', 'FontSize', 12);
% grid on;
% legend('Location', 'northeast', 'FontSize', 10);
% set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Griglie minori per maggiore dettaglio
% 
% % --- Subplot 2: Fase (in gradi) ---
% subplot(2,1,2);
% 
% % L'offset '+pi' nella fase SDOF e il segno negativo in MDOF dipendono dai tuoi dati
% % Se le curve non si allineano bene, prova a rimuovere/modificare questi offset
% plot(f_cutted, unwrap(-angle(H_pred_column_cutted)+pi), 'DisplayName', 'MDOF', 'LineWidth', 1.5, 'Color', 'r');
% hold on;
% % plot(f_cutted, unwrap(angle(frf_approx_cutted)), 'DisplayName', 'SDOF', 'LineWidth', 1.5, 'Color', 'r');
% plot(f_cutted, unwrap(angle(H_exp_cutted)), 'DisplayName', 'Experimental', 'LineWidth', 1.5, 'Color', 'b');
% 
% xlabel('Frequenza [Hz]', 'FontSize', 11);
% ylabel('Fase [°]', 'FontSize', 11);
% title('Phase FRF Comparison: Models vs. Experimental', 'FontSize', 12);
% grid on;
% legend('Location', 'northeast', 'FontSize', 10);
% set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Griglie minori
% 
% % --- FINE CODICE PLOT MIGLIORATO ---


%%

% initial parameters
Rl_overall_real = 0;
Rl_overall_imag = 0;
Rh_overall_real = 0;
Rh_overall_imag = 0;

x0 = [Rl_overall_real, Rl_overall_imag, Rh_overall_real, Rh_overall_imag];
Resid_pre = x0;

frf_SDOF = [frf_approx_cutted(1:end/2); frf_formula(end/2:end-1)];

% Optimization options
options = optimoptions('lsqnonlin', ...
'Display', 'iter-detailed', ... % Mostra più dettagli ad ogni iterazione (utile per debugging)
'MaxIterations', 2000, ...     % Aumenta leggermente le iterazioni se non converge
'MaxFunctionEvaluations', 4000, ...  % Aumenta le valutazioni della funzione
'FunctionTolerance', 1e-10, ... % Rendi più stringente
'StepTolerance', 1e-15, ...    % Rendi più stringente
'OptimalityTolerance', 1e-10); ... % Rendi più stringente

% recalling optimization func
x_optim = lsqnonlin(@(x) error_func_overall(x, H_exp_cutted, H_pred_column_cutted, omega_cutted), x0, [], [], options);

% saving optimized parameters
Resid_post = x_optim;

Rl_overall_real = 0;
Rl_overall_imag = 0;
Rh_overall_real = -0.4;
Rh_overall_imag = 0.01;

% recomputing formula
H_pred_column_cutted_optimized = H_pred_column_cutted + (((Rl_overall_real + 1i*Rl_overall_imag)./(omega_cutted).^2) + (Rh_overall_real + 1i*Rh_overall_imag)).';

%%

% --- INIZIO CODICE PLOT MIGLIORATO ---

figure('Name', 'Confronto FRF Modulo e Fase - Modelli', 'NumberTitle', 'off'); % Titolo generale della finestra

% --- Subplot 1: Modulo (scala logaritmica) ---
subplot(2,1,1);
semilogy(f_cutted, abs(H_exp_cutted), 'DisplayName', 'Experimental', 'LineWidth', 1.5, 'Color', 'b'); % Verde MATLAB
hold on;
% semilogy(f_cutted, abs(frf_formula), 'DisplayName', 'SDOF', 'LineWidth', 1.5, 'Color', 'r'); % Arancio MATLAB
% semilogy(f_cutted, abs(frf_prony_cutted), 'DisplayName', 'MDOF Prony', 'LineWidth', 1.5, 'Color', 'c'); % Blu MATLAB
semilogy(f_cutted, abs(H_pred_column_cutted_optimized), 'DisplayName', 'MDOF Ibrahim', 'LineWidth', 1.5, 'Color', 'k'); % Blu MATLAB

xlabel('Frequenza [Hz]', 'FontSize', 11);
ylabel('Modulo |H(f)|', 'FontSize', 11);
title('Module FRF Comparison: MDOF Ibrahim vs. Experimental', 'FontSize', 20);
grid on;
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Griglie minori per maggiore dettaglio

% --- Subplot 2: Fase (in gradi) ---
subplot(2,1,2);

% L'offset '+pi' nella fase SDOF e il segno negativo in MDOF dipendono dai tuoi dati
% Se le curve non si allineano bene, prova a rimuovere/modificare questi offset
plot(f_cutted, H_exp_cutted_phase, 'DisplayName', 'Experimental', 'LineWidth', 1.5, 'Color', 'b');
hold on;
% plot(f_cutted, (angle(frf_SDOF)), 'DisplayName', 'SDOF', 'LineWidth', 1.5, 'Color', 'r');
% plot(f_cutted, unwrap(-angle(frf_prony_cutted)), 'DisplayName', 'MDOF Prony', 'LineWidth', 1.5, 'Color', 'c');
plot(f_cutted, unwrap(-angle(H_pred_column_cutted_optimized)+pi), 'DisplayName', 'MDOF Ibrahim', 'LineWidth', 1.5, 'Color', 'k');

xlabel('Frequenza [Hz]', 'FontSize', 11);
ylabel('Fase [rad]', 'FontSize', 11);
title('Phase FRF Comparison: MDOF Ibrahim vs. Experimental', 'FontSize', 20);
grid on;
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Griglie minori

% --- FINE CODICE PLOT MIGLIORATO ---