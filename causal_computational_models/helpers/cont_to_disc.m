function [ disc_signal ] = cont_to_disc( cont_signal, num_bins )
%Converts a continuous signal to a discrete signal according to num bins
%you want in the discrete signal.
maxi = max(max(cont_signal));
mini = min(min(cont_signal));
least_count = (maxi-mini)/num_bins;
disc_signal = least_count*floor(cont_signal/least_count);
end

