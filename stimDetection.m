function [stim_out] = stimDetection(stim_in)
% This functions detects the start and end times of stimulations 
% stim = double with raw data from a stimulation measurement 
% int_stim = 2x(amount of stimulations) matrix, start and end times of the detected stimulations 
stim_out = stim_in ; 

for i = 1:size(stim_in,2)
    stim = stim_in{2,i} ; 
    if ~isempty(stim)
    
    absol = abs(stim) ;         % absolute value 
    mov = movstd(absol,10000) ;   % standard deviation of moving standard deviation
    m = std(mov) ; 
        
    [~,locs] = findpeaks(mov, ...
        'MinPeakProminence', (mean(mov))*3);  % detect start and ends 
    
    % Remove zeros, they mess with the ranges 
    if isempty(locs) == 0 
    [zks,~] = find(~mov(1:locs(1))) ;
    n = length(zks) ; 
    div = (n*(n+1))/2 ; 
    if zks ~= 0 
        if isequal(sum(zks), div) 
            mov = mov(locs(1):end) ;
        end
    end 
        % end of measurement not correct 
    [fks,~] = find(~mov(locs(end):end)) ;
    f = length(mov) ; 
    g = length(mov)-length(fks); 
    div = (f*(f+1))/2 - (g*(g+1))/2 ; 
    if fks ~= 0 
        if isequal(sum(fks), div) 
            mov = mov(1:locs(end)) ;
        end 
    end 
        % re-calculate peaks >> should not find peak as half is removed 
    [~,locs] = findpeaks(mov, ...
        'MinPeakProminence',(mean(mov)*3));
        % re-calculate range 
    m = std(mov) ; 
    end 
   
    % STIM OR NOSTIM
    if m < 10 % boundary for nostim range 
    
        if isempty(locs) == 1 
            startstim = [] ; 
            endstim = [] ; 
        elseif isempty(locs) == 0 
            startstim = [] ; 
            endstim = [] ;
        else 
            startstim = NaN ; 
            endstim = NaN ; 
        end 
        int_stim = {[startstim; endstim]} ; % put inside cell to move 
    
    elseif m > 10 % stim present >> standard deviation high enough for stimulation
    
        if isempty(locs) == 1 % cont stim 
            startstim = 1 ;    % beginning stim is one 
            endstim = length(mov) ;  % end stim is end 
    
        elseif isempty(locs) == 0
            if (rem(length(locs),2) ~= 0)==1
                % Start: low std left and high std right  
                % Start: high std left and low std right  
                startstim = zeros(((length(locs)+1)/2),1) ; 
                endstim = zeros(((length(locs)+1)/2),1) ;
    
                for l = 1:length(locs)
                    % Take range of 1000 for window 
                    left = abs(stim(locs(l)-1000:locs(l))) ; 
                    right = abs(stim(locs(l):locs(l)+1000)) ; 
                    for ll = 1:(length(locs)+1)/2
                        if sum(left) < sum(right) %start 
                            startstim(ll) = locs(l) ; 
                        elseif sum(right) < sum(left) %end 
                            endstim(ll) = locs(l) ;
                        else 
                            startstim = NaN ; 
                            endstim = NaN ; 
                        end 
                    end 
                end 
            
            % find if start or end had zero >> cannot be both (odd)
            % missing start is order in ascending order, zero front 
            % missing end is replace zero with last indice 
    
                if sum(endstim==0)~=0 
                    endstim(:,end) = length(mov) ; 
                elseif sum(startstim==0)~=0 
                    sort(startstim) ; 
                else 
                    startstim = NaN ; 
                    endstim = NaN ; 
                end 
            % even locs > stim started and ended 
            elseif (rem(length(locs),2) == 0) ==1
                startstim = zeros((length(locs)/2),1) ; 
                endstim = zeros((length(locs)/2),1) ;
    
                for l = 1:length(locs)
                    % Take range of 1000 for window  
                    left = abs(stim(locs(l)-1000:locs(l))) ; 
                    right = abs(stim(locs(l):locs(l)+1000)) ; 
                    for ll = 1:length(locs)/2
                        if sum(left) < sum(right) %start 
                            startstim(ll) = locs(l) ; 
                        elseif sum(right) < sum(left) %end 
                            endstim(ll) = locs(l) ; 
                        else 
                            startstim = NaN ; 
                            endstim = NaN ; 
                        end 
                    end 
                end 
            else 
                startstim = NaN ; 
                endstim = NaN ; 
            end 
        else 
            startstim = NaN ; 
            endstim = NaN ; 
        end 
        int_stim = {[startstim endstim]} ; % put inside cell to move 
    else 
        int_stim = {NaN(1,2)}; 
    end 
    stim_out(4,i) = int_stim ;
    else 
    stim_out(4,1) = {NaN(1,2)}; 
    end 
end 
end 


