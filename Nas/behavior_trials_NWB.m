

nwb_file = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41\YutaMouse41.nwb';
nwb2 = nwbRead(nwb_file);



behavior_perfect = load('20170505_396um_0um_merge.track.behavior.mat');
events_perfect = behavior_perfect.behavior.events;

events_perfect.trials          = % 1x221 cell   : structs
events_perfect.map             = % 1x10 cell    : structs
events_perfect.trialConditions = % 1x221 double : 
events_perfect.trialIntervals  = % 221x2 double :          [trials.start_time.data.load trials.stop_time.data.load]
events_perfect.conditionType   = % 1x10 cell    : strings (central, wheel)                                trials.vectordata.get('condition').data



trials = nwb2.intervals_trials;


trials.vectordata.get('both_visit').data.load
trials.vectordata.get('condition').data
trials.vectordata.get('error_run').data.load
trials.vectordata.get('stim_run').data.load




