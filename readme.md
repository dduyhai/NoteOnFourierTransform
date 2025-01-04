$\newcommand\lr[1]{\left(#1\right)}$
$\newcommand\Lr[1]{\left[#1\right]}$
$\newcommand\LR[1]{\left\{#1\right\}}$
$\newcommand\abs[1]{\left\vert #1 \right\vert}$
$\newcommand\Exp[1]{\exp\lr{#1}}$
$\newcommand\Floor[1]{\left\lfloor#1\right\rfloor}$

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
  H\lr{f} & = \int_{-\infty}^{+\infty} h\lr{t} \Exp{-2\pi \jmath f t} dt,  \\
  h\lr{t} & = \int_{-\infty}^{+\infty} H\lr{f} \Exp{2\pi \jmath f t} df.
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
Denote $f_c := \dfrac{1}{2 \Delta t}$ Nyquist critical frequency.
Sampling theorem says that if a continuous function $h\lr{t}$, sampled at an interval $\Delta t$,
happens to be __bandwidth limited__ to frequencies strictly smaller in magnitude than $f_c$
(aka $H\lr{f} = 0$ for all $\abs{f} \geqslant f_c$), then the function $h\lr{t}$ is completely
determined by its samples $h_n$.

In fact, $h\lr{t}$ is given explicitly by the formula
$$
  h\lr{t} = \Delta t \sum_{n = -\infty}^{+\infty} h_n
    \dfrac{\sin\lr{2\pi f_c\lr{t_c - n \Delta t}}}{\pi \lr{t - n \Delta t}}
$$

### Discrete Fourier transform

Suppose that we have $N$ values sampled consecutively at sampling frequency $F_s = \dfrac{1}{\Delta t}$:
$$
  h_k = h\lr{t_k = k \Delta t}, \qquad k = 0, \dots, N - 1.
$$
To make think simpler, we assume that $N = F_s$, aka., we sample only one period of sampling process.
According to sampling theorem, if $F_s$ is big enough, one can reasonably estimate $H\lr{f}$ for
all $\abs{f} < \dfrac{1}{2} F_s$. Therefore, it is plausible to seek estimations of Fourier transform 
$H\lr{f}$ only at the discrete values
$$
  f_m := \dfrac{m}{N \Delta t} = \dfrac{m}{N} F_s, \qquad m = -M, \dots, M, 
  \qquad M := \Floor{\dfrac{N}{2}}.
$$
Then $H\lr{f_m}$ is approximated by
$$
\begin{align*}
  H\lr{f_m} 
  & = \int_{-\infty}^{+\infty} h\lr{t} \Exp{-2 \pi \jmath f_m t} dt \\
  & \approx \Delta t \sum_{k = 0}^{N - 1} h_k \Exp{-2 \pi \jmath t_k f_m}, 
    \qquad m = -M, \dots, M.
\end{align*}
$$

Let us denote __discrete Fourier transform__ of $N$ points $\LR{h_n}$ as
$$
  H_n := \sum_{k = 0}^{N - 1} h_k \Exp{-2 \pi \jmath \dfrac{k n}{N}}, \quad n = 0, \dots, N - 1.
$$
We will see how array $\LR{H_n, n = 0, \dots, N - 1}$ approximate $\LR{H\lr{f_m}, m = -M, \dots, M}$.
From the formula of $H_n$, we see that 

* For $n = 0, \dots, M$: $\colorbox{yellow}{$\Delta t H_n \approx H\lr{f_n}$}$.
* $H_{n}$ is periodic in $n$ with period $N$ aka. $\colorbox{orange}{$H_{-n} = H_{N - n}$}$.
  It means that $\colorbox{yellow}{$H\lr{f_{-m}} \approx \Delta t H_{N - m}$}$ for $m = 1, \dots, M$.

With such observation, we have:
* For $m = 0, \dots, M$: $H\lr{f_m} \approx \Delta t H_m$.
* For $m = -M, \dots, -1$: $H\lr{f_m} \approx \Delta t H_{N + m} =: \Delta t H_n$ for $n = N - M, \dots, N - 1$.
* For $n = 0, \dots, M$: $H_n \approx F_s H\lr{f_n}$.
* For $n = M + 1, \dots, N - 1$: $H_n \approx F_s H\lr{f_{-N + n}} = F_s H\lr{-f_{N - n}}$ 
  ($ = F_s H\lr{f_{N - n}} \approx H_{N - n}$ in case of real $h\lr{t}$).

