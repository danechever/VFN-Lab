PM = daq.createSession('ni');
addAnalogInputChannel(PM,'Dev1',0,'Voltage');

PM.rate = 10;
PM.DurationInSeconds = 1/10;
