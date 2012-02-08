function [tfStackLow, tfStackHigh] = SplitChannels(tfStack)

% - If 'tfStack' is memory mapped, extract it

tfStackLow = double(bitand(tfStack, 255));   % Lowest 8 bits
tfStackHigh = double(bitshift(tfStack, -8)); % Highest 8 bits

