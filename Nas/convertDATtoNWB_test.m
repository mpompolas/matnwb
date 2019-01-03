




% Create an example nwb with the m120_200secs file



date = datetime(2018, 3, 1, 12, 0, 0);
session_start_time = datetime(date,'Format','yyyy-MM-dd''T''HH:mm:SSZZ',...
    'TimeZone','local');
nwb = nwbfile( 'source', 'Example File from m120_200secs', ...
    'session_description', 'a test NWB File', ...
    'identifier', 'm120_200secs', ...
    'session_start_time', session_start_time);


% Add channel info
device_labels = {};

for i =1:96
    device_labels{i} = 'PFC';
end
for i =97:192
    device_labels{i} = 'TEO';
end





udevice_labels = unique(device_labels, 'stable');

variables = {'id', 'x', 'y', 'z', 'imp', 'location', 'filtering', ...
    'group', 'group_name'};
for i_device = 1:length(udevice_labels)
    device_label = udevice_labels{i_device};
    
    nwb.general_devices.set(device_label,...
        types.core.Device());
    
    nwb.general_extracellular_ephys.set(device_label,...
        types.core.ElectrodeGroup(...
        'description', 'a test ElectrodeGroup', ...
        'location', 'unknown', ...
        'device', types.untyped.SoftLink(['/general/devices/' device_label])));
    
    ov = types.untyped.ObjectView(['/general/extracellular_ephys/' device_label]);
    
    elec_nums = find(strcmp(device_labels, device_label));
    for i_elec = 1:length(elec_nums)
        elec_num = elec_nums(i_elec);
        if i_device == 1 && i_elec == 1
            tbl = table(int64(1), NaN, NaN, NaN, NaN, {'PFC'}, {'filtering'}...
                , ov, {'electrode_group'},'VariableNames', variables);
        else
            tbl = [tbl; {int64(elec_num), NaN, NaN, NaN, NaN,...
                'PFC', 'filtering', ov, 'electrode_group'}];
        end
    end        
end



% add the |DynamicTable| object to the NWB file using the name |'electrodes'| (not flexible)

tbl.Properties.Description = 'my description';
electrode_table = util.table2nwb(tbl);
nwb.general_extracellular_ephys.set('electrodes', electrode_table);


ov = types.untyped.ObjectView('/general/extracellular_ephys/electrodes');

electrode_table_region = types.core.DynamicTableRegion('table', ov, ...
    'description', 'all electrodes',...
    'data', [1 height(tbl)]');



%% Add data
electrical_series = types.core.ElectricalSeries(...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', 30000., ... % Hz
    'data', raw_data.F,...
    'electrodes', electrode_table_region,...
    'data_unit','V');

nwb.acquisition.set('ECoG', electrical_series);



% electrical_series = types.core.ElectricalSeries(...
%     'timestamps', (1:1000)/200, ...
%     'starting_time_rate', 200., ... % Hz
%     'data', randn(10, 1000),...
%     'electrodes', electrode_table_region,...
%     'data_unit','V');


nwbExport(nwb, 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\raw_data_no_header_m120_200secs.nwb')





















