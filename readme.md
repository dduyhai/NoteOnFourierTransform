$\newcommand\lr[1]{\left(#1\right)}$
$\newcommand\Lr[1]{\left[#1\right]}$
$\newcommand\LR[1]{\left\{#1\right\}}$
$\newcommand\abs[1]{\left\vert #1 \right\vert}$
$\newcommand\Exp[1]{\exp\lr{#1}}$

# Using `fft` correctly

In this project, I try to briefly given summary on Fourier transform and how to correctly use MATLAB
`fft` function which calculates discrete Fourier transform (DFT) of a vector.

## Fourier transform
Some physical processes can be represented as a function of time $h\lr{t}$ (aka. time domain)
or as function of frequency $H\lr{f}$ (aka. frequency domain).
Here $t$ is time measured in second (or s for short) and $f$ is measured in Hz.
The two functions $h$ and $H$ can be converted back and forth using Fourier transforms:
$$
\begin{align}
  H\lr{f} & = \int_{-\infty}^{+\infty} h\lr{t} \Exp{2\pi \jmath f t} dt,  \\
  h\lr{t} & = \int_{-\infty}^{+\infty} H\lr{f} \Exp{-2\pi \jmath f t} df.
\end{align}
$$
Such conversion is denoted as $h\lr{t} \Leftrightarrow H\lr{f}$ and some of its properties are

* time scaling:
$
  h\lr{a t} \Leftrightarrow \dfrac{1}{\abs{a}} H\lr{\dfrac{f}{a}},
$
* frequency scaling:
$
  H\lr{b f} \Leftrightarrow \dfrac{1}{\abs{b}} h\lr{\dfrac{t}{b}}, 
$
* time shifting:
$
  h\lr{t - t_0} \Leftrightarrow H\lr{f} \Exp{2 \pi \jmath f t_0},
$
* frequency shifting:
$
  H\lr{f - f_0} \Leftrightarrow h\lr{t} \Exp{-2 \pi \jmath f_0 t}.
$

The Parseval's theorem says that the total power in a signal is the same regardless computing domain:
$$
  \int_{-\infty}^{+\infty} \abs{h\lr{t}}^2 dt = \int_{-\infty}^{+\infty} \abs{H\lr{f}}^2 df.
$$ 

The one-sided power spectral density of the function $h$ as
$$
  P_h\lr{f} = \abs{H\lr{f}}^2 + \abs{H\lr{-f}}^2, \qquad
  0 \leqslant f < +\infty.
$$

## Fourier transform of discretely sampled data

In most common situation, time-domain function $h\lr{t}$ is sampled (aka. its value is recorded) at
evenly spaced intervals in time.
Let $\Delta t$ the time interval between consecutive samples and the sequence of sampled values is
$$
  h_n := h\lr{t_n := n \Delta t}, \qquad n = -\infty, \dots, -1, 0, 1, \dots, +\infty.
$$ 
The time interval $\Delta t$ is called the sampling rate.  

### Sampling theorem
Denote $f_c := \dfrac{1}{\Delta t}$ Nyquist critical frequency.
Sampling theorem says that if a continuous function $h\lr{t}$, sampled at an interval $\Delta t$,
happens to be __bandwidth limited__ to frequencies strictly smaller in magnitude than $f_c$
(aka $H\lr{f} = 0$ for all $\abs(f) \geqslant f_c$), then the function $h\lr{t}$ is completely
determined by its samples $h_n$.

In fact, $h\lr{t}$ is given explicitly by the formula
$$
  h\lr{t} = \Delta t \sum_{n = -\infty}^{+\infty} h_n
    \dfrac{\sin\lr{2\pi f_c\lr{t_c - n \Delta t}}}{\pi \lr{t - n \Delta t}}
$$

### Discrete Fourier transform

Suppose that we have $N$ consecutive sampled values 
$$
  h_k = h\lr{t_k = k \Delta t}, \qquad k = 0, \dots, N - 1.
$$
Let us seek estimations of Fourier transform $H\lr{f}$ only at the discrete values
$$
  f_n = \dfrac{n}{N \Delta t} =: \dfrac{n}{N} F_s, \qquad n = -\dfrac{N}{2}, \dots, \dfrac{N}{2}.
$$
Here $F_s := \dfrac{1}{\Delta t}$ is called sampling frequency.
Then $H\lr{f_n}$ is approximated by
$$
  H\lr{f_n} = \int_{-\infty}^{+\infty} h\lr{t} \Exp{2 \pi \jmath f_n t} dt
  \approx \Delta t \sum_{k = 0}^{N - 1} h_k \Exp{2 \pi \jmath t_k f_n}
$$
Let us denote $H_n := \sum_{k = 0}^{N - 1} h_k \Exp{2 \pi \jmath \dfrac{k n}{N}}$
and call it as the discrete Fourier transformation of the $N$ points $\LR{h_k}$.
From the formula of $H_n$, we have $H_{n}$ is periodic in $n$ with period $N$,
$H_{-n} = H_{N - n}$ and hence $H_{-\frac{N}{2}} = H_{\frac{N}{2}}$.
With such properties we can shift $n$ from range $\Lr{-\dfrac{N}{2}, \dfrac{N}{2}}$
to range $\Lr{0, N - 1}$ with the __convention__:

* zero frequency corresponds to $n = 0$,
* positive frequencies $0 < f < f_c$ correspond to $1 \leqslant n \leqslant \dfrac{N}{2} - 1$,
* negative frequencies $-f_c < f < 0$ correspond to $\dfrac{N}{2} + 1 \leqslant n \leqslant N - 1$,
* the value $n = \dfrac{N}{2}$ corresponds to _both_ frequencies $-f_c$ and $f_c$.

Similarly, the formula for the discrete inverse Fourier transform which recovers
the set of $\LR{h_k}$ exactly from the $\LR{H_n}$ is
$$
  h_k = \dfrac{1}{N} \sum_{n = 0}^{N - 1} H_n \Exp{-2 \pi \jmath \dfrac{k n}{N}}.
$$
Hence, the discrete form of Parseval's theorem is
$$
  \sum_{k = 0}^{N - 1} \abs{h_k}^2 = \dfrac{1}{N} \sum_{n = 0}^{N - 1} \abs{H_n}^2.
$$
