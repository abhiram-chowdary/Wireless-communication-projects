clear all;
fname = ['city.p2m'];
% Read the file
fid = fopen(fname, 'r');
data = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
data = data{1}(6:end);
data = strjoin(data, '\n');
data = strsplit(data, '\n');
Lines = data;
num = str2double(Lines{1});
Lines = Lines(2:end);
li = Lines;
i = str2double(strsplit(Lines{1}, ' '));
cir = {};
while i(1) < num
    a = {};
    for j = 1:i(2)
        val = str2double(strsplit(Lines{j+1}, ' '));
        a{end+1} = val;
    end
    cir{end+1} = a;
    Lines = Lines(i(2)+2:end);
    i = str2double(strsplit(Lines{1}, ' '));
    if isempty(Lines)
        break;
    end
end
if i(1) == num
    a = {};
    for j = 1:i(2)
        val = str2double(strsplit(Lines{j+1}, ' '));
        a{end+1} = val;
    end
    cir{end+1} = a;
    
    % Process the data
    % ...
    
    disp(length(cir));
    cir = cir(~cellfun('isempty', cir));
    disp(length(cir));
end
channel_impulse_response = {};
for cluster = 1:length(cir)
    chr = [];
    for path = 1:length(cir{cluster})
        phase = cir{cluster}{path}(2);
        toa = cir{cluster}{path}(3);
        power = cir{cluster}{path}(4);
        chr(end+1) = complex(10^(power/10) * cos(phase), power * sin(phase)) * exp(-1i * 2 * pi * toa);
    end
    channel_impulse_response{end+1} = chr;
end
cir = channel_impulse_response;
%%
cir
