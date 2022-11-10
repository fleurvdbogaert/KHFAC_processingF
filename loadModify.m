function [stim,pres,name,data] = loadModify(fs_p,fs_s)
% Data is loaded into cells and modified. 
% Only one file can be selected at a time 
% Deleted dir outp

%% Variables
au = 10 ;   % Conversion factor AU to cmH2O
f_filt = 18;    % [Hz] low-pass filter frequency for pressure 

%% Select input file and load data 
[name,path] = uigetfile('*.*', ...
    'Select ONE FILE to put in') ;
name = {name} ; 

% Display selected file in command window
if isequal(name{1,1},0)
   error('File selection canceled') ;
else 
   disp('File selected: '), disp(name{1,1}) ; 
   data = importdata(fullfile(path,name{1,1})) ;
end

%% Select output directory 
if ~exist('Calculations', 'dir')
   mkdir 'Calculations'
end

%% Modify data 
% Stimulation 
num_stim = size(data.stimulation,2) ; % amount of stimulation measurements
if num_stim ==0 
    error('Your stimulation data is empty. Check your data files.')
end 
stim = cell(2,num_stim) ;
for n = 1:num_stim
    stim(1,n) = {data.stimulation(:,n)} ; 
    stim(2,n) = stim(1,n) ; %no adjustments to data 
    stim(3,n) = {fs_p} ; 

    mod = stim{2,n} ; 

    t = (0:numel(mod)-1)/fs_s ;

    % Plot
    figure 
    plot(t,mod, '-', 'LineWidth', 1, 'Color','#80B3FF'); 
    title(strcat('Stimulation channel ',num2str(n))); 
    xlabel('Time [s]', 'FontSize', 10);
    ylabel('Stimulation Voltage [AU]','FontSize', 10);
    hold on 

    % Ask what to display
    answer = questdlg(strcat('Do you want to display stimulation channel ' ...
        , num2str(n),' ?'), ...
        'Select stimulation channels', ...
        'Yes','No','No');
    switch answer
        case 'Yes'
            disp(strcat('stimulation channel ', num2str(n),' will be displayed.'))
        case 'No'
            stim(:,n) = {[]} ; 
            disp(strcat('stimulation channel ', num2str(n),' will NOT be displayed.'))
    end   
    close 
end 
       
% Pressure 
num_pres = size(data.pressure,2) ; % amount of pressure measurements 
if num_pres ==0                     % error message, empty files 
    error('Your pressure data is empty. Check your data files.')
end 
pres = cell(2,num_pres) ; % first raw, second modified 
for nn = 1:num_pres
    pres(1,nn) = {data.pressure(:,nn)} ; 
    pres(3,nn) = {fs_s} ; 

    % Calibration
    pres(2,nn) = {((pres{1,nn}+31.4)/1.13)/au} ; % -31.4: adjust if necessarry 

    % Notch filter 50Hz
    net_filter = designfilt('bandstopiir','FilterOrder',2, ...  % create notch filter for 50 Hz (48-52 Hz)
                       'HalfPowerFrequency1',49.9,'HalfPowerFrequency2', ...
                       50.1,'DesignMethod','butter','SampleRate',fs_p);  
    pres(2,nn) = {filtfilt(net_filter,pres{2,nn})} ; 

    % Lowpass filter
    [b,a]=butter(4,f_filt/(0.5*fs_p));  
    pres(2,nn) = {filtfilt(b,a,pres{2,nn})};

    mod = pres{2,nn} ; 
    t = (0:numel(mod)-1)/fs_p ;
    % Plot
    figure 
    plot(t,mod, '-', 'LineWidth', 1, 'Color','#80B3FF');
    title(strcat('pressure channel ',num2str(nn))); 
    xlabel('Time [s]', 'FontSize', 10);
    ylabel('Pressure [cmH_2O]','FontSize', 10);
    hold on 

    % Ask what to display
    answer = questdlg(strcat('Do you want to display pressure channel ' ...
        , num2str(nn),' ?'), ...
        'Select pressure channels', ...
        'Yes','No','No');
    switch answer
        case 'Yes'
            disp(strcat('pressure channel ', num2str(nn),' will be displayed.'))
        case 'No'
            pres(:,nn) = [] ;  % if is is the last cell line, they will ve removed 
            disp(strcat('pressure channel ', num2str(nn),' will NOT be displayed.'))
    end 
    close 
end 

