


%% Just the part that needs to be embedded



load('C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge\right_before_7.mat','trials')

x_coord = 8;
y_coord = 10;
z_coord = 9;

bins = 200;

disp('READY TO ENTER: 7')




%% Load the behavior file

my_behavior = load('C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge\20170505_396um_0um_merge.track.behavior.mat');
my_behavior = my_behavior.behavior;

bins = 200;


%% normalize positions to template
c=1;

uniqueConditions = unique((my_behavior.events.trialConditions));

for iCondition = 1:length(uniqueConditions)
    
    map{iCondition}=[];
    t_conc=[];
    
    
    iTrialsInCondition = find(my_behavior.events.trialConditions == uniqueConditions(iCondition));

    
    for iTrial = 1:length(iTrialsInCondition)
        
        selected_trial = my_behavior.events.trials{iTrialsInCondition(iTrial)}; % I set this for faster typing
        
        % Add a check here if there are x,y,z
        t_conc = [selected_trial.timestamps, selected_trial.x, selected_trial.y, selected_trial.z, 20*(selected_trial.timestamps - selected_trial.timestamps(1))]; % Check what this 20 is
%       t_conc = [trials{iCondition}{iTrial}(:,:),                    20*(trials{iCondition}{iTrial}(:,1)-trials{iCondition}{iTrial}(1,1))];

        if length(t_conc)>=bins
            while length(t_conc)>bins+1
                di = pdist(t_conc);
                s = squareform(di);
                s(find(eye(size(s))))=nan;
                [a b] = min(s(:));
                [coords blah] = find(s==a);
                t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
                t_conc(coords(2),:) = [];
            end
            t_conc_all(iTrial,:,:) = t_conc;
        end
    end
    if length(iTrialsInCondition)>0
        map{iCondition} = squeeze(median(t_conc_all(:,:,:),1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign to the right behavior field
    my_behavior.events.map{iCondition}.x = map{iCondition}(:,2);
    my_behavior.events.map{iCondition}.y = map{iCondition}(:,3);
    my_behavior.events.map{iCondition}.z = map{iCondition}(:,4);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear t_conc_all
    
%     disp('finding mapping...')
%     for iTrial =1:length(trials{iCondition})  % all trial types (rotations)
%         for iPoint = 1:length(trials{iCondition}{iTrial})
%             [a b] = min(nansum(abs([trials{iCondition}{iTrial}(iPoint,1)      -map{iCondition}(:,1),      ...   % Timestamp
%                                     trials{iCondition}{iTrial}(iPoint,x_coord)-map{iCondition}(:,x_coord),...   % X_COORDINATE
%                                     trials{iCondition}{iTrial}(iPoint,z_coord)-map{iCondition}(:,z_coord),...   % Z_COORDINATE 
%                                     trials{iCondition}{iTrial}(iPoint,y_coord)-map{iCondition}(:,y_coord),...   % Y_COORDINATE
%                                    (trials{iCondition}{iTrial}(iPoint,1)-trials{iCondition}{iTrial}(1,1))*50-map{iCondition}(:,1),...  % penalty for time differences
%                  40*(iPoint./length(trials{iCondition}{iTrial})*length(map{iCondition}) - (1:length(map{iCondition})))'])'));     % penalty for order differences
%             mapping{iCondition}{iTrial}(iPoint,:) = [map{iCondition}(b,1:end) b trials{iCondition}{iTrial}(iPoint,1)];
%         end
%     end
    
    disp('finding mapping...')
    for iTrial =1:length(iTrialsInCondition)  % all trial types (rotations)
        selected_trial = my_behavior.events.trials{iTrialsInCondition(iTrial)}; % I set this for faster typing
        for iPoint = 1:length(selected_trial.timestamps)
            [a b] = min(nansum(abs([selected_trial.timestamps(iPoint)-map{iCondition}(:,1),...   % Timestamp
                                    selected_trial.x(iPoint)-map{iCondition}(:,2),...   % X_COORDINATE
                                    selected_trial.y(iPoint)-map{iCondition}(:,3),...   % Y_COORDINATE 
                                    selected_trial.z(iPoint)-map{iCondition}(:,4),...   % Z_COORDINATE
                                   (selected_trial.timestamps(iPoint)-selected_trial.timestamps(1))*50-map{iCondition}(:,1),...  % penalty for time differences
                                    (40*(iPoint./length(selected_trial.timestamps)*length(map{iCondition}) - (1:length(map{iCondition})))')])'));     % penalty for order differences
            mapping{iCondition}{iTrial}(iPoint,:) = [map{iCondition}(b,1:end) b selected_trial.timestamps(iPoint)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Assign to the right behavior field
        my_behavior.events.trials{iTrialsInCondition(iTrial)}.mapping = mapping{iCondition}{iTrial}(:,6);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    
end

disp('READY TO ENTER: 8')


%% incorporate this below to reformat things and double check trial assignments...
% load behav
% 
c=1;
clear m mm tt
for iCondition=1:length(uniqueConditions)
    iTrialsInCondition = find(my_behavior.events.trialConditions == uniqueConditions(iCondition));
    if ~isempty(iTrialsInCondition)
        iCondition{c}=trials{iCondition};
        m{c}=map{iCondition};
        mm{c}=mapping{iCondition};

        iCondition{c}=iCondition{c}(~cellfun('isempty',iCondition{c}));
        mm{c}=mm{c}(~cellfun('isempty',mm{c}));
        c=1+c;
    end
end 
dbmap=m;
mapping=mm;
trials= iCondition;

disp('READY TO ENTER: 9')






























