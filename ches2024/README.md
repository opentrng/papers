# Introduction

This directory contains the relevant files to generate results and figures for our CHES 2024 paper [Impact of the Flicker Noise on the Ring Oscillator-based TRNGs](https://tches.iacr.org/index.php/TCHES/article/view/11450), _Licinius Benea, Mikael Carmona, Viktor Fischer, Florian Pebay-Peyroula and Romain Wacquez_.

```
@article{FlickerCHES2024,
  author={Benea, Licinius and Carmona, Mikael and Fischer, Viktor and Pebay-Peyroula, Florian and Wacquez, Romain},
  title={Impact of the Flicker Noise on the Ring Oscillator-based TRNGs},
  journal={IACR Transactions on Cryptographic Hardware and Embedded Systems},
  volume={2024},
  number={2},
  year={2024},
  month={Mar.},
  pages={870–889},
  doi={10.46586/tches.v2024.i2.870-889}
}
```

# Instructions

Here are step-by-step instructions for generating the figures.

```
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib colorednoise
python generate_figures.py
```
The `generate_figures.py` file is supposed to generate directly the Figures presented in the article, whith the exception of the histograms for which the raw data cannot be made public.

# Detailed description

## Code sections

The different parts of the script are organised as follows:
1. **Lines 26 to 143**: Emulation of ASIC results 
2. **Lines 151 to 239**: Emulation of FPGA results (raw data from `fpga_counter.txt` file)
3. **Lines 246 to 421**: Elementary RO-TRNG bit generation using the emulator and entropy rate calculation for different configurations of noise (thermal and flicker)
4. **Lines 431 to 594**: Study of the autocorrelation:
   1. **Lines 442 to 467**: Study and trace the autocorrelation of a raw time data series
   2. **Lines 471 to 594**: Study of the autocorrelation function of the bit series for different configurations of noise (thermal and flicker)
5. **Lines 598 to 628**: Trace the absolute jitter for different configurations of flicker noise (1x and 10x standard quantity from ASIC)

## Main functions

1. `LSNE(x,y,a)` implementation of the Least squares normalized Error Regression algorithm from [GB06] to obtain the coefficients of the quadratic behaviour of jitter. Exemple of use: `poly = LSNE(x_values, y_values, [2, 1, 0])`  with `[2,1,0]` for the coefficients of order 2, 1 and 0 respectively
2. `ERO_bits(T1, T2, Ath, Afl, N, size)` generates bits corresponding to an Elementary RO-TRNG. Warning: the proportionality coefficients for thermal and flicker noise need to be defined or calculated first! (generally: factor_th=2, factor_fl=0.135)
variable definition:
   1. `T1/T2` periods of RO1, RO2, from which RO2 with jitter RO1 perfect - transfer of all jitter onto RO2 
   2. `Ath, Afl` thermal and flicker noise amplitudes **without their proportionality coefficients**
   3. `N` frequency divider factor
   4. `size` number of desired bits to be generated
Example of use: `ERO_bits(2e-9, 2e-9, Ath, Afl, 1000, int(1e7))` for 10M  bits generated from two 2ns identical RO Elementary RO-TRNG with standard quantities of thermal and flicker noise, frequancy divider at 1000 .
3. `H_Baudet(T1, T2, a1, N)` determination of the minimum boundary of entropy according to [BLMT11] variable definition:
   1. `T1/T2 = periods` of RO1, RO2
   2. `a1` the linear coefficient of the prabolic curve of jitter corresponding to the amplitude of thermal noise 
   3. `N` frequency divider factor
4. `entropy(bit_series, order)` determination of the entropy rate variable definition:
   1. `bit_series` input bit series
   2. `order` size of the block of bits (8 for 2^8 values)
5. `entropy_16bits(bit_series, order=16)` idem as `entropy(bit_series, order)`, but with order blocked at 16 
6. `autocorrelation(bit_series, k_all=np.linspace(0, 100, 101).astype(int))` general computation of the autocorrelation function variable definition:
   1. `bit_series` input bit series
   2. `k_all` numpy array of desired lag computation 
