function my_behavior = testing_trials

%% Load the actual behavior
my_behavior = load('C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge\20170505_396um_0um_merge.track.behavior.mat');
my_behavior = my_behavior.behavior;


behavior    = load('C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge\20170505_396um_0um_merge.track.behavior.mat');
behavior    = behavior.behavior;


%% Check the mapping fields from the location data on each trial

trials = my_behavior.events.trials;

for iTrial = 1:length(trials)
    
    % Assign the mapping
    nBins =201;
    templateVector = 1:nBins;
    target = length(trials{iTrial}.x);
    n_ent = length(templateVector);

    template = round(interp1( 1:n_ent, templateVector, linspace(1, n_ent, target) ))';
    template = circshift(template,1);
    
    template(1)   = template(2);
    template(end) = template(end-1);
    
    
%     figure(1);plot(behavior.events.trials{iTrial}.mapping)
%     hold on
%     plot(template)
%     legend ('Buzsaki Code','My Code')
%     title(['Trial: ' num2str(iTrial)])
%     grid on
%     hold off
%     w = waitforbuttonpress

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    my_behavior.events.trials{iTrial}.mapping = template;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Check the map field
nBins = 201;
templateVector = 1:nBins;

position  = zeros(3,1); % Just to get the position dimensions
map = behavior.events.map; % Original values
new_map = cell(1,length(unique(my_behavior.events.trialConditions))); % New values
uniqueConditions = unique((my_behavior.events.trialConditions));

for iCondition = 1:length(uniqueConditions)
    
    iTrialsInCondition = find(my_behavior.events.trialConditions == iCondition);

    % Initialize template vectors
    new_map{iCondition}.x = zeros(length(templateVector),1); % 201x1
    x = cell(length(templateVector),1);

    if size(position,1)>1
        new_map{iCondition}.y = zeros(length(templateVector),1); % 201x1
        y = cell(length(templateVector),1);
    end
    if size(position,1)>2
        new_map{iCondition}.z = zeros(length(templateVector),1); % 201x1
        z = cell(length(templateVector),1);
    end
    
    % Concatenate all coordinates that are assigned to each template
    % bin
    for iTrial = 1:length(iTrialsInCondition)
        for iCoordinate = 1:length(templateVector) % 1:201
            x{iCoordinate} = [x{iCoordinate} trials{1,iTrialsInCondition(iTrial)}.x(trials{1,iTrialsInCondition(iTrial)}.mapping == iCoordinate)'];
            if size(position,1)>1
                y{iCoordinate} = [y{iCoordinate} trials{1,iTrialsInCondition(iTrial)}.y(trials{1,iTrialsInCondition(iTrial)}.mapping == iCoordinate)'];
            end
            if size(position,1)>2
                z{iCoordinate} = [z{iCoordinate} trials{1,iTrialsInCondition(iTrial)}.z(trials{1,iTrialsInCondition(iTrial)}.mapping == iCoordinate)'];
            end
        end
    end

    % Take the median (?) for each bin
    for iCoordinate = 1:length(templateVector) % 1:201
        new_map{iCondition}.x(iCoordinate) = median(x{iCoordinate});
        if size(position,1)>1
            new_map{iCondition}.y(iCoordinate) = median(y{iCoordinate});
        end
        if size(position,1)>2
            new_map{iCondition}.z(iCoordinate) = median(z{iCoordinate});
        end
    end
    
%     % Plot the difference between the two computations
%    
%     figure(iCondition);
%     plot3(map{iCondition}.x,map{iCondition}.y,map{iCondition}.z);
%     hold on
%     plot3(new_map{iCondition}.x,new_map{iCondition}.y,new_map{iCondition}.z);
%     legend ('Buzsaki Code','My Code')
%     grid on
%     hold off

%     % Plot point by point
%     figure(iCondition);
%     for iPoint = 1:length(map{iCondition}.x)
%         plot3(map{iCondition}.x(iPoint),map{iCondition}.y(iPoint),map{iCondition}.z(iPoint),'*r');
%         hold on
%         plot3(new_map{iCondition}.x(iPoint),new_map{iCondition}.y(iPoint),new_map{iCondition}.z(iPoint),'og');
%         drawnow
%         legend ('Buzsaki Code','My Code')
%         grid on
%     end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_behavior.events.map = new_map;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end






