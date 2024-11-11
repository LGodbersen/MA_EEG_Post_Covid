function trl = mytrialfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

trl = [];
    % add this to the trl definition
    begsample     = event(1).sample;
    endsample     = event(2).sample;
    offset        = 0; 
trl = [round(begsample), round(endsample), offset];


