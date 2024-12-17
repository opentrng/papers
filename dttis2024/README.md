# Introduction

This directory contains the relevant files to generate results and figures for our DTTIS 2024 paper [OpenTRNG: an open-source initiative for ring-oscillator based TRNGs](https://dttis2024.org/), _Florian Pebay-Peyroula, Licinius Benea, Mikael Carmona and Romain Wacquez_.

```
@inproceedings{OpenTrngDTTIS2024,
  author={Pebay-Peyroula, Florian and Benea, Licinius-Pompiliu and Carmona, Mikael and Wacquez, Romain},
  booktitle={2024 IEEE International Conference on Design, Test and Technology of Integrated Systems (DTTIS)}, 
  title={OpenTRNG: an open-source initiative for ring-oscillator based TRNGs}, 
  year={2024},
  pages={1-6},
  keywords={Ring oscillators;Emulation;Collaboration;Throughput;Hardware;Entropy;Reliability;Security;Random number generation;Field programmable gate arrays;TRNG;PTRNG;FPGA;ASIC;Ring-oscillator},
  doi={10.1109/DTTIS62212.2024.10780212}
}
```

# Instructions

Here are step-by-step instructions for generating the figures from the paper.

1. Clone [OpenTRNG repository](https://github.com/opentrng/ptrng)
2. Set PATH to OpenTRNG in `run.sh` and execute the script to compile and generate the data on the Digilent Arty7 plateform
3. Set PATH to OpenTRNG and to generated data in the Jupyter Notebook `figures.ipynb`, then run it to draw the figures
