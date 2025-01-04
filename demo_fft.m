clear("all");

amplitudes = [0.5, 0.75, 1.0];
freqs = [0, 50, 100];
n_sampled_values = 200;
sampling_freq = 500;
fig = figure(Visible = "off");
tab_gr = uitabgroup(fig);

sampling_interval = 1.0 / sampling_freq;
sampled_times = (0 : (n_sampled_values - 1))' * sampling_interval;

for id = 1:2
    if id == 1
        sampled_values = real(exp(2 * j * pi * freqs .* sampled_times)) * amplitudes';
    else
        sampled_values = real(exp(2 * j * pi * freqs .* sampled_times)) * amplitudes';
        sampled_values = sampled_values -1.0  + 2.0* rand(size(sampled_values));
    end
    a_half = fix(n_sampled_values / 2);

    dft_coefs = fft(sampled_values);
    dft_freqs = (0 : (n_sampled_values - 1))' / n_sampled_values / sampling_interval;
    two_sided_dft_coefs = circshift(dft_coefs, a_half);
    one_sided_dft_coefs = dft_coefs(1 : (a_half + 1));
    one_sided_dft_coefs(2 : (end - 1)) = 2 * one_sided_dft_coefs(2 : (end - 1));
    one_sided_spectrum = one_sided_dft_coefs / n_sampled_values;
    if mod(n_sampled_values, 2)
        two_sided_freqs = (-a_half : 1 : a_half)' / n_sampled_values / sampling_interval;
    else
        two_sided_freqs = ((-a_half + 1) : 1 : a_half)' / n_sampled_values / sampling_interval;
    end
    one_sided_freqs = (0 : a_half)' / n_sampled_values / sampling_interval;

    if id == 1
        tab1 = uitab(tab_gr, Title = "Signal without noise");
    else
        tab1 = uitab(tab_gr, Title = "Signal with noise");
    end
    tlo = tiledlayout(tab1, 2, 3);
    ax = nexttile(tlo, 1);
    plot(ax, sampled_times, sampled_values);
    ax.Title.String = "Original signal";
    ax.XLabel.String = "Time [s]";
    ax.YLabel.String = "Signal amplitude";

    ax = nexttile(tlo, 2);
    stem(ax, (0 : (n_sampled_values - 1))', abs(dft_coefs));
    ax.Title.String = "DFT coefficients";
    ax.XLabel.String = "Index";
    ax.YLabel.String = "Coef. magnitude";

    ax = nexttile(tlo, 3);
    stem(ax, two_sided_freqs, abs(two_sided_dft_coefs)); 
    ax.Title.String = "Two-sided DFT coefficients";
    ax.XLabel.String = "Frequency [Hz]";
    ax.YLabel.String = "Coef. magnitude";

    ax = nexttile(tlo, 4);
    stem(ax, one_sided_freqs, abs(one_sided_dft_coefs));
    ax.Title.String = "One-sided DFT coefficients";
    ax.XLabel.String = "Frequency [Hz]";
    ax.YLabel.String = "Coef. magnitude"; 

    ax = nexttile(tlo, 5);
    stem(ax, one_sided_freqs, abs(one_sided_spectrum));
    ax.Title.String = "Spectrum of frequency";
    ax.XLabel.String = "Frequency [Hz]";
    ax.YLabel.String = "Spectrum magnitude"; 

end
fig.Visible = "on";
