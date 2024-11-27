%% Clear and close everything

clear all;
fclose('all');

% Check if serial port object or any communication interface object exists
serialobj = instrfind;
if ~isempty(serialobj)
    delete(serialobj)
end

clc;
clear all;
close all;

%% ---------- Serial port setting ----------
s = serialport("COM4",115200);
wsize = 512; % Number of sampling points, i.e., number of data points to acquire
fs = 240; % Sampling rate, check the setting in Arduino 

%% ---------- Sampling setting ----------
L1 = 4;
w1 = ones(1,L1)/L1;
L2 = 21;
w2 = zeros(1,L2); w2(ceil(L2/2)) = 1;
w2 = w2 - ones(1,L2)/L2;

w = conv(w1,w2);
L3 = length(w);
L4 = 32;
%freqz(w)


%%
%time_axis = 1:wsize; % Time axis of the display buffer
disp_data = nan(wsize, 1);
heart_rate_est = zeros(100,1);

% Initialize figure object
fig = figure('Name','raw data');
time_axis = single((1:wsize)'/fs);

h_plot = plot(nan,nan);
hold on;
t_plot = plot(nan, nan, 'o');

% process global defination
proc_window = zeros(1, L3);
diff_window = zeros(1, 2);
smooth_window = zeros(1, L4);
proc_idx = 1;
diff_idx = 1;
smth_idx = 1;

% main procedure
i = 1;
m = 1;
while (1)
    % terminal condition
    k=get(gcf,'CurrentCharacter');
    if k=='p'
        filename = "pic.jpg";
        saveas(fig, filename);
        break;
    end
    if k=='q', break; end


    line = readline(s);
    raw_data = str2double(line);

    [flt_data, proc_window, proc_idx] = data_process(raw_data, w, proc_window, proc_idx);

    % Add data to display buffer
    if i < wsize
        disp_data(i) = flt_data;
    else
        disp_data = [disp_data(2:end); flt_data]; % first in first out
    end

    % Update figure plot
    if(mod(i, 10) == 0)
        [max_value, max_idx] = findpeaks(disp_data);
        for j = 1:length(max_idx)
            if max_value(j) < 50
                max_value(j) = nan;
                max_idx(j) = nan;
            end
        end
        max_value = rmmissing(max_value);
        max_idx = rmmissing(max_idx);
        
        dif = [max_idx;0] - [0;max_idx];
        dif = dif(2:length(dif)-1);
        rate = sum(dif) / length(dif);

        set(h_plot, 'xdata', time_axis, 'ydata', disp_data);
        set(t_plot, 'xdata', max_idx/fs, 'ydata', max_value);
        
        %rate = 1;
        str = strcat('the heart rate: ', num2str(fs / rate * 60, 3), 'BPM');
        title(str);
        xlabel('Time');
        ylabel('Quantized value');
        drawnow;
        
        if i > wsize && mod(i, 30) == 0 && m <= 100 
            heart_rate_est(m) = rate;
            disp(m);
            m = m+1;
        end
    end
    
    if i == 1 tic
    end
    if i == 2400 toc 
    end
    i = i+1;
end

%%
% Disconnect the serial port object from the serial port
clear s;



%%
function [flt_data, process_window, idx] = data_process(raw_data, w, process_window, idx)
    L = length(w);
    process_window(mod(idx, L)+1) = raw_data;
    d_rt = 0;
    for i = 1:L
        d_rt = d_rt + process_window(mod(idx+i-1, L)+1) * w(i);
    end
    flt_data = d_rt;
    idx = idx+1;
end




