# Introduction

This directory contains the relevant files to generate results and figures for our DTTIS 2024 paper [OpenTRNG: an open-source initiative for ring-oscillator based TRNGs](https://dttis2024.org/), _Florian Pebay-Peyroula, Licinius Benea, Mikael Carmona and Romain Wacquez_.

```
@article{DTTISOpenTRNG2024,
  author={Pebay-Peyroula, Florian and Benea, Licinius and Carmona, Mikael and Wacquez, Romain},
  title={OpenTRNG: an open-source initiative for ring-oscillator based TRNGs},
  journal={IEEE International conference on Design Test and Technology of Integrated Systems},
  year={2024},
  month={Oct.}
}
```

# Instructions

Here are step-by-step instructions for generating the figures from the paper.

1. Clone [OpenTRNG repository](https://github.com/opentrng/ptrng)
2. Set PATH to OpenTRNG in `run.sh` and execute the script to compile and generate the data on the Digilent Arty7 plateform
3. Set PATH to OpenTRNG and to generated data in the Jupyter Notebook `figures.ipynb`, then run it to draw the figures
