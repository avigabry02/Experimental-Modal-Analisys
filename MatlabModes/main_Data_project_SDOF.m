clc
close all
clear all 


%% Loading of the files
% selecting main folder where are all the data subfolder sets
% mainFolder = 'C:\Users\Utente\Desktop\LAUREA_MAGISTRALE\SEMESTER_2\DATA ANALYSIS FOR MECHANICAL SYSTEM IDENTIFICATION\PROJECT\experimental_data';
% mainFolder = 'C:\Users\Carlo Tundo\Desktop\PoliMi\Master 1st year\DATA ANALYSIS\assignement\modal analysis\test modal analysis';
mainFolder = '/Users/gabri/Documents/università /PoliMi/Lessons/First Year/Second Semester/Data Analysis for Mechanical Systems Identification/Project/MatlabModes/test modal analysis';

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
% % optional
% set = 1; % default
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
H1_C = cell(nset, 1);
H1_complex = cell(nset, 1);
Coherence = cell(nset, 1);
neg = 0; % to cancel out negative frequencies

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

omega = f.*2*pi;
omega(1) = 0.2;
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

%% -----------  SDOF  ------------ 
%% circle calculation for SDOF

for ii = 1:nset 
    % H1-COMPLEX
    H1_complex{ii}.out1 = Cspectra{ii}.acc_Cs1./Pspectra{ii}.input_Ps;
    H1_complex{ii}.out2 = Cspectra{ii}.acc_Cs2./Pspectra{ii}.input_Ps;
    H1_complex{ii}.out3 = Cspectra{ii}.acc_Cs3./Pspectra{ii}.input_Ps;
    H1_complex{ii}.out4 = Cspectra{ii}.acc_Cs4./Pspectra{ii}.input_Ps;
    H1_complex{ii}.out5 = Cspectra{ii}.acc_Cs5./Pspectra{ii}.input_Ps;
end

% selecting frf to plot circles
npeaks = 6;
ii = idx_set;
jj = idx_acc;
cell_h1_complex = {H1_complex{ii}.out1, H1_complex{ii}.out2, H1_complex{ii}.out3, H1_complex{ii}.out4, H1_complex{ii}.out5};

% select frf
frf = cell_h1_complex{jj};
real_part_frf = real(frf);
imag_part_frf = imag(frf);

% pay attention to this!!! see with code lines after this if he is finding
% every peak correctly 
[pks, idx_pks] = findpeaks(abs(frf), "MinPeakHeight", 25, 'MinPeakDistance', 400);

figure(10)
semilogy(f2, abs(frf))
hold on 
scatter(f2(idx_pks), pks)
grid on 

% select the first 6 modes
pks = pks(1:6);
idx_pks = idx_pks(1:6);

% define the ranges
freq0 = f2(idx_pks);
min_f = freq0 - 20;
max_f = freq0 + 20;

figure
for zz = 1:npeaks
    % finding index 
    [~, idx_min] = min(abs(f2 - min_f(zz)));
    [~, idx_max] = min(abs(f2 - max_f(zz)));
    
    % select for each mode its range
    frf_ranged = frf(idx_min:idx_max);
    real_ranged_frf = real_part_frf(idx_min:idx_max);
    imag_ranged_frf = imag_part_frf(idx_min:idx_max);
    

    len = 1:length(frf_ranged);
    len_fine = linspace(1, length(frf_ranged), 5000);

    Re_spline = spline(len, real(frf_ranged), len_fine);
    Im_spline = spline(len, imag(frf_ranged), len_fine);

    subplot(2,3,zz)
    plot(Re_spline, Im_spline, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Spline interpolante');
    hold on
    scatter(real_ranged_frf, imag_ranged_frf, 'r');
    title(sprintf('Mode %d', zz));
    xlabel('Re')
    ylabel('Im')
    axis equal
    grid on
end


%% computing sdof modal analysis 
max_eig = 6; % we don't have more than 6 eigenfreq

FRF_first = cell(nset, 1);
COHERENCE_first = cell(nset, 1);
freq_range = cell(nset, 1);
f0 = cell(nset, 1);
w0 = cell(nset, 1);
csi0 = cell(nset, 1);
f_M = cell(nset, 1);
omega_M = cell(nset, 1);
FRF = cell(nset, 1);
good_sensors = cell(nset, 1);
G_num = cell(nset, 1);
Initial_param = cell(nset, 1);
Final_param = cell(nset, 1);

vector_h1 = cell(nset, 1);
vector_h1_complex = cell(nset, 1);
vector_phase = cell(nset, 1);
vector_coh = cell(nset, 1);


for ii = 1:nset
    vector_h1{ii} = {H1{ii}.abs_out1; H1{ii}.abs_out2; H1{ii}.abs_out3; H1{ii}.abs_out4; H1{ii}.abs_out5};
    vector_h1_complex{ii} = {H1_complex{ii}.out1; H1_complex{ii}.out2; H1_complex{ii}.out3; H1_complex{ii}.out4; H1_complex{ii}.out5};
    vector_phase{ii} = {Phase{ii}.phase_acc1; Phase{ii}.phase_acc2; Phase{ii}.phase_acc3; Phase{ii}.phase_acc4; Phase{ii}.phase_acc5};
    vector_coh{ii} = {Coherence{ii}.gamma_acc1; Coherence{ii}.gamma_acc2; Coherence{ii}.gamma_acc3; Coherence{ii}.gamma_acc4; Coherence{ii}.gamma_acc5};

    % Division in ranges
    FRF_first{ii} = cell(max_eig, 1);
    COHERENCE_first{ii} = cell(1, max_eig);
    f_M{ii} = cell(max_eig, 1);
    omega_M{ii} = cell(max_eig, 1);
    FRF{ii} = cell(max_eig, 1);

    for jj = 1:n_acc
    % everything by eye - tutto a occhio!!!
        [~, idx1] = min(abs(f2 - 237));
        [~, idx2] = min(abs(f2 - 257));
        FRF_first{ii}{1}(:, jj) = vector_h1_complex{ii}{jj}(idx1:idx2);
        COHERENCE_first{ii}{1}(:, jj) = vector_coh{ii}{jj}(idx1:idx2);
        freq_range{ii}(:, 1) = f2(idx1:idx2);
    
    
        [~, idx3] = min(abs(f2 - 318));
        [~, idx4] = min(abs(f2 - 338));
        FRF_first{ii}{2}(:, jj) = vector_h1_complex{ii}{jj}(idx3:idx4);
        COHERENCE_first{ii}{2}(:, jj) = vector_coh{ii}{jj}(idx3:idx4);
        freq_range{ii}(:, 2) = f2(idx3:idx4);
    
        [~, idx5] = min(abs(f2 - 382));
        [~, idx6] = min(abs(f2 - 402));
        FRF_first{ii}{3}(:, jj) = vector_h1_complex{ii}{jj}(idx5:idx6);
        COHERENCE_first{ii}{3}(:, jj) = vector_coh{ii}{jj}(idx5:idx6);
        freq_range{ii}(:, 3) = f2(idx5:idx6);
    
        [~, idx7] = min(abs(f2 - 586));
        [~, idx8] = min(abs(f2 - 606));
        FRF_first{ii}{4}(:, jj) = vector_h1_complex{ii}{jj}(idx7:idx8);
        COHERENCE_first{ii}{4}(:, jj) = vector_coh{ii}{jj}(idx7:idx8);
        freq_range{ii}(:, 4) = f2(idx7:idx8);
    
        [~, idx9] = min(abs(f2 - 610));
        [~, idx10] = min(abs(f2 - 630));
        FRF_first{ii}{5}(:, jj) = vector_h1_complex{ii}{jj}(idx9:idx10);
        COHERENCE_first{ii}{5}(:, jj) = vector_coh{ii}{jj}(idx9:idx10);
        freq_range{ii}(:, 5) = f2(idx9:idx10);
    
        [~, idx11] = min(abs(f2 - 796));
        [~, idx12] = min(abs(f2 - 816));
        FRF_first{ii}{6}(:, jj) = vector_h1_complex{ii}{jj}(idx11:idx12);
        COHERENCE_first{ii}{6}(:, jj) = vector_coh{ii}{jj}(idx11:idx12);
        freq_range{ii}(:, 6) = f2(idx11:idx12);
        % 
        % [~, idx13] = min(abs(f2 - 945));
        % [~, idx14] = min(abs(f2 - 965));
        % FRF_first{ii}{7}(:, jj) = vector_h1_complex{ii}{jj}(idx13:idx14);
        % COHERENCE_first{ii}{7}(:, jj) = vector_coh{ii}{jj}(idx13:idx14);
        % freq_range{ii}(:, 7) = f2(idx13:idx14);
        % 
        % [~, idx15] = min(abs(f2 - 1013));
        % [~, idx16] = min(abs(f2 - 1033));
        % FRF_first{ii}{8}(:, jj) = vector_h1_complex{ii}{jj}(idx15:idx16);
        % COHERENCE_first{ii}{8}(:, jj) = vector_coh{ii}{jj}(idx15:idx16);
        % freq_range{ii}(:, 8) = f2(idx15:idx16);
        % 
        % [~, idx17] = min(abs(f2 - 1191));
        % [~, idx18] = min(abs(f2 - 1211));
        % FRF_first{ii}{9}(:, jj) = vector_h1_complex{ii}{jj}(idx17:idx18);
        % COHERENCE_first{ii}{9}(:, jj) = vector_coh{ii}{jj}(idx17:idx18);
        % freq_range{ii}(:, 9) = f2(idx17:idx18);
    end

    %% Finding eigenfrequencies
    f0_sum = cell(max_eig, 1);
    f0{ii} = zeros(max_eig, 1);
    good_sensors{ii} = cell(max_eig, 1);
    
    for ll = 1:max_eig      % for each mode
        k = 0;        % counter for the sensor in which the i-th mode has good coherence
        temp = []; 
        for jj = 1:n_acc   % for each sensor
            
            % selecting peak value and its index
            [max_value, idx_max] = max(abs(FRF_first{ii}{ll}(:, jj)));
    
            if (COHERENCE_first{ii}{ll}(idx_max, jj) >= 0.8)
                % when coherence is higher than 0.85
                freq_peak = freq_range{ii}(idx_max, ll);
                temp = [temp; freq_peak];
                k = k + 1;
                good_sensors{ii}{ll}(k) = jj; % saving sensor index for each mode
            end
        end
        % average of the frequencies
        f0_sum{ll} = temp;
        f0{ii}(ll) = mean(f0_sum{ll}); 
       
    end
    % omegas at which we have resonances
    w0{ii} = f0{ii}.*2.*pi;
    
    
    %% Better and symmetric division around w0
    f_min = f0{ii} - 3;
    f_max = f0{ii} + 3;
    omega_min = 2*pi*f_min;
    omega_max = 2*pi*f_max;
    
    % and the corrisponding indexes are
    idx_min = [];
    idx_max = [];
    idx_f0 = [];
    
    % f_M{ii}{ll} = cell(max_eig, 1);
    % omega_M{ii}{ll} = cell(max_eig, 1);
    % FRF{ii}{ll} = cell(max_eig, 1);
    
    for ll = 1:max_eig
        [~, idx_min(ll)] = min(abs(f2 - f_min(ll)));
        [~, idx_max(ll)] = min(abs(f2 - f_max(ll)));
        [~, idx_f0(ll)] = min(abs(f2 - f0{ii}(ll)));
    
        f_M{ii}{ll} = f2(idx_min(ll) : idx_max(ll));
        omega_M{ii}{ll} = f_M{ii}{ll}.*2*pi;
    
        for jj = 1:n_acc
            FRF{ii}{ll}(:, jj) = vector_h1_complex{ii}{jj}(idx_min(ll) : idx_max(ll));
        end
    end
    
    
    %% Finding damping ratios 
    csi0_sum = cell(max_eig, 1);
    csi0{ii} = zeros(max_eig, 1);
    
    for ll = 1:max_eig
        temp = []; 
        for jj = good_sensors{ii}{ll}
            [~, idx_peak] = min(abs(f_M{ii}{ll} - f0{ii}(ll)));
            half_power_value = abs(FRF{ii}{ll}(idx_peak, jj)) / sqrt(2);
            idx = find(abs(FRF{ii}{ll}(:, jj)) >= half_power_value);
            idx_left = idx(1);
            idx_right = idx(end);
            f1_damping = f_M{ii}{ll}(idx_left);
            f2_damping = f_M{ii}{ll}(idx_right);
            
            damping_val = ((f2_damping*2*pi)^2 - (f1_damping*2*pi)^2) / (4 * (f0{ii}(ll)*2*pi)^2);
            temp = [temp; damping_val]; 
    
        end
        csi0_sum{ll} = temp;
        csi0{ii}(ll) = mean(csi0_sum{ll}); 
    end
    
    
    
    %% Initial estimation of parameters: MODES' CONSTANTS (A, Rl and Rh)
    % preallocating variables
    A0 = cell(max_eig, 1);
    Rl0_real = cell(max_eig, 1);
    Rl0_imag = cell(max_eig, 1);
    Rh0_real = cell(max_eig, 1);
    Rh0_imag = cell(max_eig, 1);
    
    for ll = 1:max_eig    % for each mode
        A0{ll} = zeros(n_acc, 1);
        Rl0_real{ll} = zeros(n_acc, 1);
        Rl0_imag{ll} = zeros(n_acc, 1);
        Rh0_real{ll} = zeros(n_acc, 1);
        Rh0_imag{ll} = zeros(n_acc, 1);

        % start with the co-located
        jj = ii;
        [~, idx_peak] = min(abs(f_M{ii}{ll} - f0{ii}(ll)));
        
        phi_ii = sqrt(imag(FRF{ii}{ll}(idx_peak, jj))*2*csi0{ii}(ll)); % in terms of accelerations!!!
        A_col = phi_ii^2;
        A0{ll}(jj) = A_col;

        for jj = 1:n_acc     % for each other sensor non co-located
                
            if jj~=ii
                phi_jj = (imag(FRF{ii}{ll}(idx_peak, jj))*2*csi0{ii}(ll))/phi_ii; % in terms of accelerations!!!
                A_noncol = phi_jj*phi_ii; 

                A0{ll}(jj) = A_noncol;
            end

                
            % % initial values for A constants
            % A0{ll}(jj) = -imag(FRF{ii}{ll}(idx_peak, jj))*2*csi0{ii}(ll)*(w0{ii}(ll))^2;
            
            % initial values for residuals
            Rl0_real{ll}(jj) = 0;
            Rl0_imag{ll}(jj) = 0;
            Rh0_real{ll}(jj) = 0;
            Rh0_imag{ll}(jj) = 0;
        end
    end
    
    
    %% Optimization 
    % preallocation of frf numerical and parameters variables
    G_num{ii} = cell(max_eig, 1);
    Initial_param{ii} = cell(max_eig, 1);
    Final_param{ii} = cell(max_eig, 1);
    
    for ll = 1:max_eig
        % initial parameters
        x0 = [w0{ii}(ll); csi0{ii}(ll); A0{ll}; Rl0_real{ll}; Rl0_imag{ll}; Rh0_real{ll}; Rh0_imag{ll}];
        Initial_param{ii}{ll} = x0;
        
        % Optimization options
        options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...         % Mostra output a ogni iterazione
        'MaxIterations', 1000, ...     % Numero massimo di iterazioni
        'MaxFunctionEvaluations', 2000, ...  % Max valutazioni funzione
        'FunctionTolerance', 1e-9, ... % Tolleranza su funzione obiettivo
        'StepTolerance', 1e-12, ...    % Tolleranza sui passi
        'OptimalityTolerance', 1e-9, ... % Tolleranza sul gradiente
        'Algorithm', 'trust-region-reflective');  % Algoritmo (default)

        lb = zeros(length(x0), 1);
        ub = [2*w0{ii}(ll); 2*csi0{ii}(ll); 2*A0{ll}; 100*ones(5,1); 100*ones(5,1); 100*ones(5,1); 100*ones(5,1)];
    
        % recalling optimization func
        x_optim = lsqnonlin(@(x) error_func2(x, FRF{ii}{ll}, omega_M{ii}{ll}, n_acc), x0, [], [], options);
    
        % saving optimized parameters
        Final_param{ii}{ll} = x_optim;
        
        % now we create our numerical tf
        % r = length(A0{ii});
        for kk = 1:n_acc
            w = x_optim(1);
            csi = x_optim(2);
            A      = x_optim(3           : 2+n_acc);
            Rl_re  = x_optim(3+n_acc         : 2+2*n_acc);
            Rl_im  = x_optim(3+2*n_acc      : 2+3*n_acc);
            Rh_re  = x_optim(3+3*n_acc       : 2+4*n_acc);
            Rh_im  = x_optim(3+4*n_acc       : 2+5*n_acc);
            % G_num{3}(:, kk) = ((A0{3}(kk)./(-(omega_M{3}).^2 + 2*1i*csi0(3)*w.*(omega_M{3}) + w0(3).^2)) + ((Rl0_real{3}(kk) + 1i*Rl0_imag{3}(kk))./(omega_M{3}).^2) + (Rh0_real{3}(kk) + 1i*Rh0_imag{3}(kk))).';
            G_num{ii}{ll}(:, kk) = ((A(kk).*((omega_M{ii}{ll}).^2))./(-(omega_M{ii}{ll}).^2 + 2*1i*csi*w.*(omega_M{ii}{ll}) + w.^2)) + ((Rl_re(kk) + 1i*Rl_im(kk))./(omega_M{ii}{ll}).^2) + (Rh_re(kk) + 1i*Rh_im(kk)).';
        end
    end
    

end

%% Plots
% selceting index by input
input = inputdlg('Selecting input sensor to plot: ');
ii = str2double(input{1});

sens = inputdlg('Selecting output sensor to plot: ');
jj = str2double(sens{1});

figure
subplot(2,1,1)
semilogy(f, vector_h1{ii}{jj}(1:length(f)), 'LineWidth',2)
hold on
scatter(f_M{ii}{1}, abs(G_num{ii}{1}(:, jj)))
scatter(f_M{ii}{2}, abs(G_num{ii}{2}(:, jj)))
scatter(f_M{ii}{3}, abs(G_num{ii}{3}(:, jj)))
scatter(f_M{ii}{4}, abs(G_num{ii}{4}(:, jj)))
scatter(f_M{ii}{5}, abs(G_num{ii}{5}(:, jj)))
scatter(f_M{ii}{6}, abs(G_num{ii}{6}(:, jj)))
% scatter(f_M{ii}{7}, abs(G_num{ii}{7}(:, jj)))
% scatter(f_M{ii}{8}, abs(G_num{ii}{8}(:, jj)))
% scatter(f_M{ii}{9}, abs(G_num{ii}{9}(:, jj)))
grid on 
title(compose("Amplitude of the %d-th sensor", jj));
xlabel('f')

subplot(2,1,2)
plot(f, vector_phase{ii}{jj}(1:length(f)), 'LineWidth',2)
hold on 
scatter(f_M{ii}{1}, angle(G_num{ii}{1}(:, jj)))
scatter(f_M{ii}{2}, angle(G_num{ii}{2}(:, jj)))
scatter(f_M{ii}{3}, angle(G_num{ii}{3}(:, jj)))
scatter(f_M{ii}{4}, angle(G_num{ii}{4}(:, jj)))
scatter(f_M{ii}{5}, angle(G_num{ii}{5}(:, jj)))
scatter(f_M{ii}{6}, angle(G_num{ii}{6}(:, jj)))
% scatter(f_M{ii}{7}, angle(G_num{ii}{7}(:, jj)))
% scatter(f_M{ii}{8}, angle(G_num{ii}{8}(:, jj)))
% scatter(f_M{ii}{9}, angle(G_num{ii}{9}(:, jj)))
grid on 
title(compose('Phase of the %d-th sensor', jj))
xlabel('f')

%% recostruction of FRF 
% preallocation of variables
G_num_overall = cell(nset,1);

for ii = 1:nset
    G_num_overall{ii} = cell(max_eig, 1);
    for ll = 1:max_eig
        for kk = 1:n_acc
            w = Final_param{ii}{ll}(1);
            csi = Final_param{ii}{ll}(2);
            A      = Final_param{ii}{ll}(3           : 2+n_acc);
            Rl_re  = Final_param{ii}{ll}(3+n_acc         : 2+2*n_acc);
            Rl_im  = Final_param{ii}{ll}(3+2*n_acc      : 2+3*n_acc);
            Rh_re  = Final_param{ii}{ll}(3+3*n_acc       : 2+4*n_acc);
            Rh_im  = Final_param{ii}{ll}(3+4*n_acc       : 2+5*n_acc);
            % G_num_overall{ii}{ll}(:, kk) = ((A(kk)./(-(omega).^2 + 2*1i*csi*w.*(omega) + w.^2)) + ((Rl_re(kk) + 1i*Rl_im(kk))./(omega).^2) + (Rh_re(kk) + 1i*Rh_im(kk))).';
            % G_num_overall{ii}{ll}(:, kk) = ((A(kk)./(-(omega).^2 + 2*1i*csi*w.*(omega) + w.^2)) + (Rh_re(kk) + 1i*Rh_im(kk))).';
            G_num_overall{ii}{ll}(:, kk) = (((A(kk).*(omega).^2)./(-(omega).^2 + 2*1i*csi*w.*(omega) + w.^2)) ).';

        end
    end
end

% sum of all the modes 
frf_approx = zeros(length(f),1);
ii = idx_set;
kk = idx_acc;
for ll = 1:max_eig
    frf_approx = frf_approx + G_num_overall{ii}{ll}(:, kk);
end

% % plotting values
% figure
% subplot(2,1,1)
% semilogy(f(1:22500), vector_h1{ii}{kk}(1:22500), 'LineWidth', 1.2,'Color', 'b')
% hold on 
% semilogy(f(1:22500), abs(frf_approx(1:22500)), 'LineWidth', 1.2, 'Color', 'r')
% xlabel('frequency [Hz]')
% ylabel('abs')
% title(['FRF comparison for Input accelerometer ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
% grid on 
% legend({'FRF', 'Approx'}, 'Location', 'northwest');
% 
% subplot(2,1,2)
% plot(f(1:22500),vector_phase{ii}{kk}(1:22500), 'LineWidth', 1.2, 'Color', 'b')
% hold on
% plot(f(1:22500),angle(frf_approx(1:22500)), 'LineWidth', 1.2, 'Color', 'r')
% grid on
% xlabel('frequency [Hz]')
% ylabel('angle [rad]')
% title(['Phase Comparison for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
% legend({'FRF', 'Approx'}, 'Location', 'northwest');

%% improving SDOF with overall residuals
ii = idx_set;
jj = idx_acc;

% slicing of vectors
[~, idx_min_frf] = min(abs(f - 240));
[~, idx_max_frf] = min(abs(f - 840));

frf_cutted = frf_approx(idx_min_frf:idx_max_frf);
omega_cutted = omega(idx_min_frf:idx_max_frf);

% initial parameters
Rl_overall_real = 0;
Rl_overall_imag = 0;
Rh_overall_real = 0;
Rh_overall_imag = 0;

x0 = [Rl_overall_real, Rl_overall_imag, Rh_overall_real, Rh_overall_imag];
Resid_pre = x0;

% Optimization options
options = optimoptions('lsqnonlin', ...
'Display', 'iter', ...         % Mostra output a ogni iterazione
'MaxIterations', 1000, ...     % Numero massimo di iterazioni
'MaxFunctionEvaluations', 2000, ...  % Max valutazioni funzione
'FunctionTolerance', 1e-9, ... % Tolleranza su funzione obiettivo
'StepTolerance', 1e-12, ...    % Tolleranza sui passi
'OptimalityTolerance', 1e-9, ... % Tolleranza sul gradiente
'Algorithm', 'trust-region-reflective');  % Algoritmo (default)

% recalling optimization func
x_optim = lsqnonlin(@(x) error_func_overall(x, vector_h1_complex{ii}{jj}(idx_min_frf:idx_max_frf), frf_cutted, omega_cutted), x0, [], [], options);

% saving optimized parameters
Resid_post = x_optim;

Rl_overall_real = x_optim(1);
Rl_overall_imag = x_optim(2);
Rh_overall_real = x_optim(3);
Rh_overall_imag = x_optim(4);

% recomputing formula
frf_formula = frf_cutted + (((Rl_overall_real + 1i*Rl_overall_imag)./(omega_cutted).^2 + (Rh_overall_real + 1i*Rh_overall_imag))).';

% plotting
figure
subplot(2,1,1)
semilogy(f(idx_min_frf:idx_max_frf), vector_h1{ii}{kk}(idx_min_frf:idx_max_frf), 'LineWidth', 1.2,'Color', 'b')
hold on 
semilogy(f(idx_min_frf:idx_max_frf), abs(frf_formula), 'LineWidth', 1.2, 'Color', 'r')
xlabel('frequency [Hz]')
ylabel('abs')
title(['FRF comparison for Input accelerometer ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
grid on 
legend({'FRF', 'Approx'}, 'Location', 'northeast');

subplot(2,1,2)
plot(f(idx_min_frf:idx_max_frf),vector_phase{ii}{kk}(idx_min_frf:idx_max_frf), 'LineWidth', 1.2, 'Color', 'b')
hold on
plot(f(idx_min_frf:idx_max_frf),angle(frf_formula), 'LineWidth', 1.2, 'Color', 'r')
grid on
xlabel('frequency [Hz]')
ylabel('angle [rad]')
title(['Phase Comparison for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
legend({'FRF', 'Approx'}, 'Location', 'northeast');