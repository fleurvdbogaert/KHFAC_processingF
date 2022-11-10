function [T,Tav_pres,Ts] = calcExport(stim,pres,name,data,fs_p)

%% Avoid error and get outcomes if no stimulation channel is displayed
nostimcheck= stim(5,:);
emptyCells = cellfun(@isempty,nostimcheck);
for c=1:length(emptyCells)
if sum(emptyCells) == length(emptyCells) & isempty(pres{5,c}) == 0
    stim(1,c)= {data.stimulation(:,c)} ; 
    stim(2,c)= stim(1,c);
    stim(3,c)= {fs_p}; 
    stim(4,c)= {[0,0]};
    stim(5,c)={[0,1,length(data.stimulation)]};
end
end
%% Get outcomes 
for i = 1:size(stim,2); % amount of channels used stimulation 
    smod = stim{2,i} ; 
    sfs = stim{3,i} ; 
    time = stim{4,i} ;
    smat = stim{5,i} ; 
    %% Create plot 
if ~isempty(smat)
    figure
    hold on; 
    t = linspace(1,size(smod,1),size(smod,1)); 
    dir = strcat(pwd,'/Calculations') ; 
    file = name{1,1}; 
    filename = sprintf('%s/Stim%d%s%s.xlsx', dir,i,file,datestr(now,'yyyymmdd_HH_MM_ss')); % state filename of the output file
         
    for ii = 1:size(pres,2) % amount of channels used pressure  
        pmod = pres{2,ii} ; 
        pfs = pres{3,ii} ;
        tme = pres{4,ii} ;
        pmat = pres{5,ii} ; 

        tp = (0:numel(pmod)-1)/(pfs); % Time vector 

    if ~isempty(pmat)==1
        hold on ; 
        subplot(size(pres,2),1,ii) ; xlim([0 size(smod,1)/sfs]) ; 

        str = '#80B3FF';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

        for tt = 1:size(time) 
            patch([time(tt,:)./sfs time(tt,2)./sfs time(tt,1)./sfs] ...
                ,[-20 -20 20 20],color,'Edgecolor',color);%'Color','#80B3FF') ; 
        end 
        st1 = sprintf('START CONTRACTION CHANNEL %d',i); 
        st2 = sprintf('END CONTRACTION CHANNEL %d',i); 

        hold on 
        xline(tme(:,1)./pfs,'g-',st1) ; 
        xline(tme(:,2)./pfs,'r-',st2) ; 
        plot(tp, pmod, 'b-', 'LineWidth', 3 );%,'Color','#80B3FF'); % mod data
        xlabel('Time [s]', 'FontSize', 10);
        hold on
        set(gcf, 'Position',  [200, 200, 1000, 400])      % make a rectangular figure
 %% Individual contractions - Outcomes 
        Xp = {'no contraction', 'contraction'} ;
        Xs = {'no stimulation', 'stimulation'} ;

        % Labels of intervals. 
        label = Xp(1+pmat(:,1))' ;  % 1. labels
        
        % Start time, end time and duration of the contraction
        starttime = (pmat(:,2)-1)./pfs ; starttime(1,1) = 0;      % 2.start
        peaktime = (pmat(:,3)-1)./pfs ;       % 3.stop
        duration = peaktime-starttime ; % 4.duration
        
        % Start pressure, end pressure and absolute height of
        % contraction, all values are for one contraction. 
        startpres = pmod(round(pmat(:,2),0),1) ; 
        peakpres = pmod(round(pmat(:,3),0),1) ; 
        absheight = peakpres - startpres ; % 5.height
        
        % Slope of the contraction per contraction
        slope = absheight./duration ;    % 6.slope 
        
        % Integral calculated per contraction
        mask = pmat(:,2) <= t & t <= pmat(:,3) ; %create 7x7200000 array with masks 
        modnew = repmat(abs(pmod),1,size(startpres,1))' ; %create 7x7200000 matrix with mod 
        integral = trapz(modnew.*mask,2) ;   % 7. integral: integrate over rows 
        
        %% Averages of all contractions - outcomes 
        
        numcon = sum(pmat(:,1)) ; % total number of contractions
        avdurtot = (sum(duration.*pmat(:,1)))/numcon ; % average contraction time 
        avstartpres = (sum(startpres.*pmat(:,1)))/numcon ; % average start pressure contraction
        avpeakpres = (sum(peakpres.*pmat(:,1)))/numcon ; 
        avabsheight = (sum(absheight.*pmat(:,1)))/numcon ;
        avslope = (sum(slope.*pmat(:,1)))/numcon ;    % 6.slope 
        
        %% All stimulations - outcomes  

        % Labels of intervals. 
        labels = Xs(1+smat(:,1))' ;  % 1. labels
        
        % Start time, end time and duration of the contraction
        stim_starttime = (smat(:,2)-1)./sfs ; starttime(1,1) = 0;       % 2.start
        stim_stoptime = (smat(:,3)-1)./sfs ;        % 3.stop
        stim_duration = (stim_stoptime-stim_starttime) ; % 4.duration
        % Mean and max pressure during every (non)stimulation interval
        for g=1:length(smat(:,1));
        beginp= fix(smat(g,2));
        endp= fix(smat(g,3));
        Av_pres_stiminterval(g,:) = mean(pmod(beginp:endp));
        Max_pres_stiminterval(g,:) = max(pmod(beginp:endp));
        end

        %% Contractions within stimulation interval 
      
        Number_contr_start = zeros(size(smat,1),1) ; % size amount stimulation intervals
        Number_contr_end = zeros(size(smat,1),1) ; % size amount stimulation intervals

        for ff = 1:size(smat,1) % all stimulation intervals 
            for f = 1:size(tme,1) % contraction intervals 
                
                strt = 0; 
                stp = 0; 
        
                if (tme(f,1)>= smat(ff,2) && tme(f,1)<= smat(ff,3)) == 1 %contraction starts within stim interval
                    strt = 1 ; 
                end 
        
                if (tme(f,2)>= smat(ff,2) && tme(f,2)<= smat(ff,3)) == 1 %contraction ends within stim interval
                    stp = 1 ; 
                end 
        
                if strt == 1
                    Number_contr_start(ff,1)=(Number_contr_start(ff,1)+1); %counts contractions starts within stim interval
                end
                if stp == 1
                    Number_contr_end(ff,1)=(Number_contr_end(ff,1)+1); %counts contractions ends within stim interval
                end 
            end 
        end
    
   
        %% Tables 
            T = table(label, ...
                starttime, peaktime, duration,...
                startpres, peakpres, absheight,... 
                slope, integral) ;
            Ts = table(labels, ...
                stim_starttime, stim_stoptime, stim_duration, Av_pres_stiminterval, Max_pres_stiminterval, Number_contr_start, Number_contr_end) ;
            Tav_pres = table(numcon, ...
                avdurtot, avstartpres, avpeakpres, avabsheight,... 
                avslope); 

        rng1 = sprintf('A%d:H%d',size(pmat,1)+3,size(pmat,1)+4) ;
        rng2 = sprintf('A%d:J%d',size(pmat,1)+6,size(pmat,1)+(length(stim_duration)+8)) ;

        warning('off','MATLAB:xlswrite:AddSheet'); %optional
        writetable(T,filename,'Sheet',ii);
        writetable(Tav_pres,filename,'Sheet',ii,'Range',rng1);
        writetable(Ts,filename,'Sheet',ii,'Range',rng2);

    elseif isempty(pmat)==1
    
    writetable(cell2table({'No pressure measurement'}), ...
        filename,'Sheet',ii,'Range','A1:A2') ;  
    
    end
    end     
end 
end 







