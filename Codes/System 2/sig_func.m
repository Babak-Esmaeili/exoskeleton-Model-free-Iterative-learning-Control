function out = sig_func(signal,power)

%     out = (abs(signal).^power).*sign(signal);
    out = (abs(signal).^power).*tanh(signal/0.01);

end