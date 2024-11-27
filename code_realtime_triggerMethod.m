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


%%
%time_axis = 1:wsize; % Time axis of the display buffer
disp_data = int16(nan(wsize, 1));
trig_data = int16(nan(wsize, 1));

% Initialize figure object

fig = figure('Name','raw data');
time_axis = single((1:wsize)'/fs);

h_plot = plot(nan,nan);
hold on;
t_plot = plot(nan, nan, 'o');


% process global defination
process_window = zeros(1, L3);
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

    [flt_data, process_window, proc_idx] = data_process(raw_data, w, process_window, proc_idx);

    [t_data, diff_window, diff_idx, smooth_window, smth_idx] = ...
        data_process_(flt_data, diff_window, diff_idx, smooth_window, smth_idx);

    % Add data to display buffer
    if i < wsize
        disp_data(i) = flt_data;
        trig_data(i) = t_data;
    else
        disp_data = [disp_data(2:end); flt_data]; % first in first out
        trig_data = [trig_data(2:end); t_data];
    end

    % Update figure plot
    if(mod(i, 10) == 0)
        [rate, trig_idx] = find_trigger(trig_data, 100);
        rate = fs / rate * 60;

        trig_value = zeros(1, length(trig_idx));
        for j = 1:length(trig_idx)
            trig_value(j) = disp_data(trig_idx(j));
        end
        
        set(h_plot, 'xdata', time_axis, 'ydata', disp_data);
        set(t_plot, 'xdata', trig_idx/fs, 'ydata', trig_value);
        
        %rate = 1;
        str = strcat('the heart rate: ', num2str(rate, 3), 'BPM');
        title(str);
        xlabel('Time');
        ylabel('Quantized value');
        drawnow;
    end
    
    if i == 1 tic
    end
    if i == 4800 toc 
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
    d_rt = int16(d_rt);
    flt_data = d_rt;
    idx = idx+1;
end


function [flt_data, df_window, df_idx, sm_window, sm_idx] ...
         = data_process_(data, df_window, df_idx, sm_window, sm_idx)
    % difference filter
    df_window(mod(df_idx, 2)+1) = data;
    df_rt = df_window(mod(df_idx+1, 2)+1) - df_window(mod(df_idx, 2)+1);
    df_idx = df_idx + 1;
    
    % smooth filter
    L = 32;
    sm_window(mod(sm_idx, L)+1) = df_rt * df_rt;
    sm_rt = 0;
    for i = 1:L
        sm_rt = sm_rt + sm_window(mod(sm_idx+i, L)+1);
    end
    sm_rt = sm_rt / L;
    sm_idx = sm_idx + 1;
    
    flt_data = sm_rt;
end


function [rate, idx] = find_trigger(data, threshold)
    L = length(data);
    idx = [];
    for i = L:-1:2
        trig = xor(data(i)>threshold, data(i-1)>threshold);
        incr = data(i) > data(i-1);

        if trig && incr
            idx = [idx, i];
        end
    end
    idx_ = [0 idx] - [idx 0];
    idx_ = idx_(2:length(idx_)-1);
    rate = sum(idx_) / length(idx_);
end


