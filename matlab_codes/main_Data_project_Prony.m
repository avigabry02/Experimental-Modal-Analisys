clc
close all
clear all 


%% Loading of the files

% selecting main folder where are all the data subfolder sets
%mainFolder = '\Users\enric\OneDrive\Desktop\Data analysis for mechanical systems\test modal analysis';
%mainFolder = 'C:\Users\Utente\Desktop\LAUREA_MAGISTRALE\SEMESTER_2\DATA ANALYSIS FOR MECHANICAL SYSTEM IDENTIFICATION\PROJECT\experimental_data';
mainFolder = '/Users/gabri/Documents/università /PoliMi/Lessons/First Year/Second Semester/Data Analysis for Mechanical Systems Identification/Project/MatlabModes/test modal analysis';
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
% optional
% set = 5; % default
% 
% plotting_title = [{'Input'}, {'Acc1'}, {'Acc2'}, {'Acc3'}, {'Acc4'}, {'Acc5'}];
% for jj = 1:ntest
%     figure('Name', ['Test ' num2str(jj)], 'NumberTitle', 'off');
%     plotting_cell = {test{set}.data{jj}.Dati(:, 1), test{set}.data{jj}.Dati(:, 2), test{set}.data{jj}.Dati(:, 3), test{set}.data{jj}.Dati(:, 4), test{set}.data{jj}.Dati(:, 5), test{set}.data{jj}.Dati(:, 6)};
%     for kk = 1:(n_acc + 1)
%         subplot(n_acc+1,1,kk)
%         plot(t, plotting_cell{kk}, 'LineWidth', 1.2)
%         title(plotting_title{kk})
%         xlabel ('f [Hz]')
%         ylabel('EU')
%         grid on 
%     end
%     sgtitle(['Time domain plots: set ', num2str(set)]);
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
title(['H1, H2 and Coherence for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
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
title(['Phase for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])


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


%% Prony method

h_selected_sectioned = h(3,:,:);

h_selected = squeeze(h_selected_sectioned);

%%
N_responses = length(h_selected(:,1));

N_vector = [12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144 , 156 , 168 , 180 , 192 , 204 , 216, 228];

poles_tot = [];
A_tot = [];

for N_prony = N_vector

    R_total = [];
    B_total = [];

    m = 200*N_prony;
    
    for ii=1:N_responses
        h = h_selected(ii,:);
        R_local = zeros(m,2*N_prony);
        B_local = zeros(m,1);
        for jj=1:m
            R_local(jj,:) = h(jj : jj + 2*N_prony  -1);
            B_local(jj) = h(jj + 2*N_prony);
        end
        R_total = [R_total ; R_local];
        B_total = [B_total ; B_local];
    end
    
    
    Betas = -1.*pinv(R_total)*B_total;
    
    Betas_tot = [Betas ; 1];
    
    V = roots(flip(Betas_tot));
    
    V_matrix = zeros(2*N_prony , 2*N_prony);
    exponents = 0:1:2*N_prony -1;

    for hh = 1:2*N_prony
        row = V.^exponents(hh);
        V_matrix(hh,:) = row;
    end
    
    A_local = zeros(2*N_prony,N_responses);

    for ii = 1:N_responses
        h_vector = (R_total((ii-1)*m + 1,:)).';
        A_local(:,ii) = V_matrix\h_vector;
    end

    poles_local = log(V)./dt;
    
    poles_tot = [poles_tot ; poles_local];
    A_tot = [A_tot ; A_local];

end

A_tot = A_tot.';

%% Frequency stabilization:

freq_stab = abs(poles_tot)./(2*pi);


freq_res = 0.6; %width of the frequency range for the stabilization diagram 0.1%
edges = min(freq_stab):freq_res:820;

figure
histogram(freq_stab, 'BinEdges', edges)

freq_stable_poles = cell(6, 1); %we put here all the poles that are stable in frequency
A_stable_modes = cell(6, 1);
damping_cell = cell(6, 1);
freq_range = [247.119, 328.119, 391.719, 596.319, 620.319, 806.319];


%% Damping stabilization:


for ii = 1:length(freq_range)
    kk = 0;
    for jj = 1:length(poles_tot)
        
        current_mode = A_tot(:, jj);
        current_pole = poles_tot(jj);
        current_freq = abs(current_pole)./(2*pi);

        if(current_freq >= freq_range(ii) && current_freq <= (freq_range(ii)+freq_res))
            kk = kk + 1;
            freq_stable_poles{ii, 1}(kk) = current_pole;
            A_stable_modes{ii, 1}(:, kk) = current_mode;
            %we compile the damping cell 
            current_damping = (-1)*real(current_pole)./(current_freq*2*pi);
            damping_cell{ii, 1}(kk) = current_damping;
        end
    end
end

%%
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
A_damping_stable_modes = cell(6, 1);

for ii = 1:6
    kk = 0;
    for jj = 1:length(freq_stable_poles{ii, 1}(:))

        current_mode = A_stable_modes{ii, 1}(:, jj);
        current_freq = abs(freq_stable_poles{ii, 1}(jj))./(2*pi);
        current_damping = (-1)*real(freq_stable_poles{ii, 1}(jj))./(current_freq*2*pi);

        if(current_damping >= mostFrequentInterval_damping(ii,1) && current_damping <= mostFrequentInterval_damping(ii,2))
            kk = kk + 1;

            A_damping_stable_modes{ii, 1}(:, kk) = A_stable_modes{ii, 1}(:, jj);
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

freq_damping_stable_modes = cell(6, 1);
for ii = 1:6
    modes_ = A_damping_stable_modes{ii, 1};
    poles_ = freq_damping_stable_poles{ii, 1};
    for jj = 1:size(modes_, 2)

        current_omega = abs(poles_(jj));
        current_damping = (-1)*real(poles_(jj))./(current_omega);
        freq_damping_stable_modes{ii, 1}(:,jj) = modes_(:,jj).*2.*1i.*current_omega .*sqrt(1 - current_damping^2);
    end
end

colocated = cell(6,1);
for ii = 1:6
    colocated{ii, 1} = real(sqrt(freq_damping_stable_modes{ii, 1}(3,:)));
end

SPERIAMO = cell(6, 1);
for ii = 1:6
    for jj = 1:5
        SPERIAMO{ii, 1}(jj,:) = freq_damping_stable_modes{ii, 1}(jj, :)./colocated{ii, 1};
    end
end


freq_damping_stable_modes_real = cell(6, 1);
for ii = 1:6
    freq_damping_stable_modes_real{ii, 1}(:, :) = real(SPERIAMO{ii, 1}(:, :));
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
modes_average_quad = [];
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
    modes_average_quad(2*ii - 1, :) = mean(current_mode, 2);
    modes_average_quad(2*ii, :) = mean(current_mode, 2);

end

%%
modes_average_normalized = [];
for ii = 1:12
    max_val = max(abs(modes_average_quad(ii,:)));
    modes_average_normalized(ii, :) = modes_average_quad(ii,:)./max_val;
end

% %%
% H_colocated = 0;
% H_pred_modal = 0;
% H_pred_modal_2 = 0;
% 
% % slicing of vectors
% [~, idx_min_frf] = min(abs(f - 240));
% [~, idx_max_frf] = min(abs(f - 840));
% f_cutted = f(idx_min_frf:idx_max_frf);
% 
% for ii = 1:2:12
%     H_colocated = H_colocated + ((modes_average_quad(ii, 3).*((2*pi.*f_cutted).^2))./(-(2*pi.*f_cutted).^2 + 2*1i.*(2*pi.*f_cutted).*(2*pi*f0_average(ii)).*csi_average(ii) + (2*pi*f0_average(ii))^2));
% 
% end
% %% plotting
% 
% figure
% subplot(2, 1, 1)
% semilogy(f_cutted, abs(H_pred_modal), 'LineWidth', 1.2, 'Color', 'r');
% hold on
% semilogy(f_cutted, abs(H1_C{3, 1}.out3(idx_min_frf:idx_max_frf)), 'LineWidth', 1.2, 'Color', 'b');
% % semilogy(f_new, abs(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'k');
% legend('Mdof pinv', 'Experimental', 'Location','northwest')
% 
% subplot(2, 1, 2)
% plot(f_cutted, angle(H_pred_modal)+pi, 'LineWidth', 1.2, 'Color', 'r');
% hold on
% plot(f_cutted, angle(H1_C{3, 1}.out3(idx_min_frf:idx_max_frf)), 'LineWidth', 1.2, 'Color', 'b');
% % semilogy(f_new, abs(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'k');
% legend('Mdof pinv', 'Experimental', 'Location','northwest')
% 
% colocated_row = sqrt(modes_average_quad(5,:));
% modes_average= zeros(size(modes_average_quad));
% 

% %%
% for ii=1:2*6
%     modes_average(ii,:) = modes_average_quad(ii,:)./colocated_row;
% end
% 
% for i = 1:size(modes_average, 1)
%     for j = 1:size(modes_average, 2)
%         if isreal(modes_average(i, j))
%             modes_average(i, j) = real(modes_average(i, j));
%         else
%             modes_average(i, j) = -1*imag(modes_average(i, j));
%         end
%     end
% end
% %%
% norm_modes_average= zeros(size(modes_average_quad));
% 
% for ii=1:2*6
%     norm_factor = max(abs(modes_average(ii,:)));
%     norm_modes_average(ii,:) = modes_average(ii,:)./norm_factor;
% end
% 
% 
% 
% 

%% Modal mass:

N_modes = 6;
f_new = f(4:end);
A = zeros(length(f_new),(N_modes+2));
X = zeros((2*N_modes + 2), 1);
B = zeros(N_modes,1);

B = H1_C{3, 1}.out3(1:length(f_new)); %this is Hpq where p = sensor and q = input

for ii = 1:(N_modes)

    %in this case we use 3 because the sensor is the 3-th
    column_numerator = ones(length(f_new), 1).*modes_average_quad((2*ii)-1, 3).^2;
    column_denominator = 1i.*(2*pi.*f_new') - (ones(length(f_new), 1).*poles_average((2*ii)-1));

    A(:, ii) = column_numerator./column_denominator;

end

A(:, end-1) = ones(length(f_new), 1);

A(:, end) = 1./((f_new.*2*pi).^2);


X = pinv(A) * B;  % pseudoinversa più robusta

X_2 = A\B;



%% Finding m_modal:

Q_vector = X(1:6,1);
Q_vector_2 = X_2(1:6, 1);
R_right = X(7,1);
R_left=X(8,1);
modal_masses = [];
modal_masses_2 = [];

for ii = 1 : N_modes

    modal_masses(ii) = abs((-1*(2*pi.*f0_average(ii*2-1)).^2)/(Q_vector(ii)*2*1i*f0_average(ii*2-1)*2*pi*sqrt(1-(csi_average(ii*2-1))^2)));
    modal_masses_2(ii) = abs((-1*(2*pi.*f0_average(ii*2-1)).^2)/(Q_vector_2(ii)*2*1i*f0_average(ii*2-1)*2*pi*sqrt(1-(csi_average(ii*2-1))^2)));

end


%% FRF estimation:

H_pred_modal = 0;
H_pred_modal_2 = 0;

% slicing of vectors
[~, idx_min_frf] = min(abs(f - 240));
[~, idx_max_frf] = min(abs(f - 840));
f_cutted = f(idx_min_frf:idx_max_frf);

for ii = 1:6
    H_pred_i = -((modes_average_quad(2*ii, 3)*modes_average_quad(2*ii, 3)*(2*pi.*f_cutted).^2)./(modal_masses(ii)*(-(2*pi.*f_cutted).^2 + 2*1i.*(2*pi.*f_cutted).*(2*pi*f0_average(2*ii)).*csi_average(2*ii) + (2*pi*f0_average(2*ii))^2)));
    H_pred_modal = H_pred_modal + H_pred_i;
    H_pred_i_2 = (modes_average_quad(2*ii, 3)*modes_average_quad(2*ii, 3)*((2*pi.*f).^2))./(modal_masses_2(ii)*(-(2*pi.*f).^2 + 2*1i.*(2*pi.*f).*(2*pi*f0_average(2*ii)).*csi_average(2*ii) + (2*pi*f0_average(2*ii))^2));
    H_pred_modal_2 = H_pred_modal + H_pred_i;
end


H_pred = A*X;



%% plotting

H_pred_modal_prony = H_pred_modal;

figure
subplot(2, 1, 1)
semilogy(f_cutted, abs(H_pred_modal), 'LineWidth', 1.2, 'Color', 'r');
hold on
semilogy(f_cutted, abs(H1_C{3, 1}.out3(idx_min_frf:idx_max_frf)), 'LineWidth', 1.2, 'Color', 'b');
% semilogy(f_new, abs(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'k');
legend('Mdof pinv', 'Experimental', 'Mdof \', 'Location','northwest')

subplot(2, 1, 2)
plot(f_cutted, angle(H_pred_modal_prony), 'LineWidth', 1.2, 'Color', 'r');
hold on
plot(f_cutted, angle(H1_C{3, 1}.out3(idx_min_frf:idx_max_frf)), 'LineWidth', 1.2, 'Color', 'b');
% semilogy(f_new, abs(H_pred_modal_2(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'k');
legend('Mdof pinv', 'Experimental', 'Mdof \', 'Location','northwest')


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
MDOF_modes = modes_average_normalized(1:2:end, :)';
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
% Esempio semplificato per la ricostruzione con i residui diretti
% Assumiamo che A_damping_stable_modes{mode_idx, 1}(sensor_idx, pole_idx) contenga il residuo per quel polo e sensore
% E che freq_damping_stable_poles{mode_idx, 1}(pole_idx) contenga il polo

% H_reconstructed_direct = zeros(size(f_cutted));
% 
% % Per ogni modo stabile identificato (6 modi)
% for ii = 1:6
%     % Prendi i poli e le ampiezze stabili per questo modo
%     current_poles = freq_damping_stable_poles{ii, 1};
%     current_amplitudes = A_damping_stable_modes{ii, 1}(3,:); % Assicurati che 3 sia l'indice corretto per il tuo sensore di interesse (output)
% 
%     % I poli sono a coppie coniugate, i residui pure
%     % Per ogni coppia di polo coniugato e il suo residuo associato
%     for jj = 1:length(current_poles)
%         pole_k = current_poles(jj);
%         amplitude_k = current_amplitudes(jj);
% 
%         % Contributo di questo polo alla FRF
%         H_reconstructed_direct = H_reconstructed_direct + amplitude_k ./ (1i*2*pi*f_cutted' - pole_k);
%     end
% end

% Poi confronta H_reconstructed_direct con H1_C{3, 1}.out3(idx_min_frf:idx_max_frf)
