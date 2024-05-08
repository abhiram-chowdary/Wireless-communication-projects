% Initialize cell arrays to store power, time, and phase
allPower = {};
allTime = {};
allPhase = {};
allcir={};
m=0;
mi=100000;
timevector=2048;
track=0;
%%
for n = 1:length(cir)
    i = cir{n};
    if size(i)[1] ~= 0
    track=track + 1;    
    power = [];
    time = [];
    phase = [];
    
    for j = i
        pow = j{1}(4);
        phas = j{1}(2);
        timed = j{1}(3);
        
        power = [power, pow];
        phase = [phase, phas];
        time = [time, timed];
    end
    time=round(time*1e10);
    if isempty(time)==0
    if m<max(time) 
        m=max(time);
    end
    if mi>min(time) 
        mi=min(time);
    end
    end
    cires=zeros([1,timevector]);
    atun=sqrt(power).*exp(1i*phase);
    cires(time)=atun;
    dcires=downsample(cires,2);
    allcir{track}=dcires;
    % Store the results for each iteration in cell arrays
    allPower{track} = power;
    allTime{track} = time;
    allPhase{track  } = phase;
    end
end
m
mi
% Access the results for each iteration using allPower, allTime, allPhase
% For example, allPower{1} contains the power values for the first iteration.
