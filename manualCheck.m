function [newdata] = manualCheck(olddata,tpe)
% Manually Check if you agree with the automatically found start and end.
% If the code did not work earlier, you can fully manually add all start
% and end points here. The new start and end points will be added to the
% third row of cells and will replace the old found points. If there is a
% partial start or end, zero or last value will be taken as start or end 
% respectively.  
% tpe = type of data, pressure or stimulation 
% olddata = data with automatically located 

data = olddata ; 

for i = 1:size(data,2)   % loop over all measurements included 

        raw = data{1,i} ; % raw data 
        mod = data{2,i} ; % modified data 
        fs = data{3,i} ; % put sample frequency inside cell with data.
        pnt = data{4,i} ; % points

if ~isempty(raw) == 1 
        total_time= (numel(mod)/fs);
        t = (0:(1/fs):(total_time-(1/fs)));
        % figure  
            figure 
            plot(t, mod, '-', 'LineWidth', 1, 'Color','#80B3FF'); % mod data 
            xlabel('Time [s]', 'FontSize', 10);
            hold on
            set(gcf, 'Position',  [200, 200, 1000, 400])      % make a rectangular figure

    if sum((sum(isnan(pnt))) > 1)==1 || isempty(pnt) == 1 % if fs the data did not go through the stim or cont detection
        for ii = 1:2 % needs to be '2', is about the start and the end 
    
               % determine name for pop-ups 
               if (isequal('stimulation',tpe)) == 1 
                    ylabel('Stimulation Voltage [AU]','FontSize', 10); 
                    
                    if ii == 1 
                        name = ' STIMULATION START' ; 
                    elseif ii == 2 
                        name = ' STIMULATION END' ; 
                    end 
                elseif (isequal('pressure',tpe)) == 1
                    ylabel('Pressure [cmH_2O]','FontSize', 10);
                    
                    if ii == 1 
                        name = ' CONTRACTION ONSET' ; 
                    elseif ii == 2 
                        name = ' CONTRACTION PEAK' ; 
                    end 
               end 

               if isempty(pnt) == 1
                   disp(['No start and end were detected. ' ...
                       'Check if this is correct'])
               elseif sum(isnan(pnt)) > 1
                   disp(['The code could not pass automatic detection. ' ...
                       'Please fill out start and end points manually'])
               end 
          
            % Manual adding  
            choicep = 'Yes' ; 
            
            while isequal(choicep,'Yes')
                choicep = questdlg(strcat('Do you want to add a ' ...
                    ,name,'?'), ...
                    strcat('Add',name), ...
                'Yes','No','No');
                % Handle response
                switch choicep
                    case 'Yes'
                        title(strcat('Select',name));
                        [xn, ~] = ginput(1); 
                        pnt(end+1,ii) = xn*fs;
                        xline(xn, '-',name);
                        disp(strcat('Added',name)) 
                    case 'No'
                        disp(strcat('No more',name,' added.'))
                end    
            end 
            hold on ; 
            if ~isempty(pnt) == 1 
                pnt_left = sort(rmmissing(nonzeros(pnt(:,1)))) ; % remove zeros and sort ascending
                pnt_right = sort(rmmissing(nonzeros(pnt(:,2)))) ;
            else                 
                pnt_left = [] ; % remove zeros and sort ascending
                pnt_right = [] ;
            end 
        end 
       
    elseif (isempty(raw) && sum(sum(isnan(pnt)))) == 0 % if measurements are included (all is filled with [] if not included) 
        for ii = 1:size(pnt,2) % needs to be '2', is about the start and the end
               % determine name for pop-ups 
               if (isequal('stimulation',tpe)) == 1 
                    ylabel('Stimulation Voltage [AU]','FontSize', 10); 
                    if ii == 1 
                        name = ' STIMULATION START' ; 
                    elseif ii == 2 
                        name = ' STIMULATION END' ; 
                    end 
                elseif (isequal('pressure',tpe)) == 1
                    ylabel('Pressure [cmH_2O]','FontSize', 10);
        
                    if ii == 1 
                        name = ' CONTRACTION ONSET' ; 
                    elseif ii == 2 
                        name = ' CONTRACTION PEAK' ; 
                    end 
               end 
                
                title(strcat('Labelling',name)); 
                nantest = isnan(pnt);  
                for j=1:numel(pnt(:,ii))
                    if nantest(j,ii) == 0 %no plotting of NaN 
                    xline((pnt(j,ii)/fs), 'LineWidth',2);   % grey lines 
                    else
                    end
                hold on; % hold state on 
                end
                
                % Keep or delete  
                for iii=1:numel(pnt(:,ii))
                    if nantest(iii,ii) == 0
                        xline((pnt(iii,ii)/fs), 'b', 'LineWidth',2)
                        answer = questdlg(strcat('Is the blue line a' ...
                        ,name,'?'), ...
                        'Keep or delete?', ...
                        'Keep','Delete', 'Leave Selection', 'Keep'); 
                        % Handle response
                        switch answer
                            case 'Keep'
                                xline((pnt(iii,ii)/fs), 'g',name, 'LineWidth', 2)
                                disp(strcat('Kept',name))
                            case 'Delete'
                                xline((pnt(iii,ii)/fs), 'r', 'LineWidth', 2)
                                pnt(iii,ii) = 0;
                                disp(strcat('Deleted',name))
                            case 'Leave Selection' 
                                close all force                            
                                return  
                        end
                    else
                    end    
                end 
        
                % Add extra  
                choicep = 'Yes' ; 
                
                while isequal(choicep,'Yes')
                    choicep = questdlg(strcat('Do you want to add another' ...
                        ,name,'?'), ...
                        'Add extra', ...
                    'Yes','No','No');
                    % Handle response
                    switch choicep
                        case 'Yes'
                            title(strcat('Select',name));
                            [xn, ~] = ginput(1); 
                            pnt(end+1,1:2) = NaN;
                            pnt(end,ii) = xn*fs;
                            xline(xn, '-',name);
                            disp(strcat('Added',name)) 
                        case 'No'
                            disp(strcat('No more',name,' added.'))
                    end    
                end 
                hold on ; 

                if ~isempty(pnt) == 1 
                    pnt_left = sort(rmmissing(nonzeros(pnt(:,1)))) ; % remove zeros and sort ascending
                    pnt_right = sort(rmmissing(nonzeros(pnt(:,2)))) ; 
                else                 
                    pnt_left = [] ; % remove zeros and sort ascending
                    pnt_right = [] ;
                end 
        end 
    end 
    
        if isempty(pnt_left) && isempty(pnt_right) == 1 
            pnt_left = zeros(1,1) ; 
            pnt_right = zeros(1,1) ;
        elseif isempty(pnt_left) == 1       % no onset, one end 
            pnt_left = zeros(1,1) ; 
        elseif isempty(pnt_right) == 1      % one onset, no end 
            pnt_right = max(t) ; 
        else
            if pnt_right(1)>pnt_left(1) %first end bigger than start
                if isequal(numel(pnt_right), numel(pnt_left)) % every onset followed by peak
                elseif numel(pnt_right) < numel(pnt_left)  % missing peak at end 
                    pnt_right = [pnt_right; max(t)] ;
                end 
    
            elseif pnt_left(1)>pnt_right(1) %first start bigger than end 
                if numel(pnt_left) < numel(pnt_right)      % missing onset at begin 
                    pnt_left = [0;pnt_left] ; 
                elseif isequal(numel(pnt_right), numel(pnt_left))  % missing end peak and begin onset
                    pnt_right = [pnt_right; max(t)] ;
                    pnt_left = [0;pnt_left] ; 
                    %else 
                end 
            %else 
            end  
        end             
        man = cat(2,pnt_left,pnt_right) ;
    close 
else 
    man = [] ; 
end 
data(4,i) = {man} ; 

%% Get intervals 

    time = man ; 
    t2 = linspace(1,size(mod,1),size(mod,1)); 
    if ~isempty(raw) % measurement channel used

        %% Get intervals
        num = size(time,1); 
        if sum(time(1,:)) == 0   % correct for fs to time mistake earlier in code 
            mat = [0 1 t2(end)] ; 
        else
            if time(1,1) == 0 
                time(1,1) = 1 ;  % zero does not exist in fs domain and cannot be used for indexing
            end 
 
        % Intervals
        if num > 1
            if time(1,1) == 1                   % start activation
                startint = time(:,2)+1 ;        % act end = int start 
                endint = vertcat(time(2:end,1)-1,t2(end)) ; % act start = int end  
                if time(end) == t2(end)          % end with act 
                    startint = time(1:end-1,2)+1 ; % act end = int start 
                    endint = time(2:end,1)-1 ;  % act start = int end 
                end 
            elseif time(end) == t2(end)          % end with act
                startint = vertcat(1, time(1:end-1,2)+1) ; % act end = int start 
                endint = time(:,1)-1 ;          % act start = int end  
            else % no partial activation
                startint = vertcat(1, time(:,2)+1) ; % act end = int start 
                endint = vertcat(time(:,1)-1, t2(end)) ; % act start = int end              
            end 
        elseif num == 1 % single %(:,x) should not be necessary really >> add error? In case it is larger?
            if time(1,1) == 1 % partial at start 
                if (time(end)==t2(end))==1 || sum(sum(time))==1 % full or none
                   startint = []; 
                   endint = [] ; 
                else 
                   startint = time(:,2) ; 
                   endint = t2(end) ; 
                end 
            elseif time(1,2) == t2(end) % partial at end 
                startint = 1 ; 
                endint = time(:,1)-1 ; 
            else % no partial
                startint = vertcat(1, time(:,2)-1) ; 
                endint = vertcat(time(:,1)+1,t2(end)) ; 
            end 
        end 
       
        act = horzcat(ones(size(time,1),1),time) ;
        int = horzcat(zeros(size(startint,1),1), startint,endint) ;
        
        ca = vertcat(act,int) ; 
        [~,idx] = sort(ca(:,2));            % sort just the first column
        mat = ca(idx,:);                    % sort the whole matrix using the sort indices
        end 
    elseif isempty(raw)
        data(5,i) = {[]};
        continue ; 
    else 
        %error message? your raw data i contains an abnormal value. Klick
        %... in the workspace, then click the first row of the ith column.
        %What do you see? 
    end 

    data(5,i) = {mat} ; 
end
newdata = data ;
end  



