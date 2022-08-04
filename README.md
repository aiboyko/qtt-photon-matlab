# qtt-photon-matlab
2D wave scattering solver on astronomically large grids using QTT decomposition


## Citation
If you find this code useful please cite

```
@INPROCEEDINGS{8262141,
  author={Boyko, Alexey I. and Oseledets, Ivan V. and Gippius, Nikolai A.},
  booktitle={2017 Progress In Electromagnetics Research Symposium - Spring (PIERS)}, 
  title={Towards solving lippmann-schwinger integral equation in 2D with polylogarithmic complexity with quantized tensor train decomposition}, 
  year={2017},
  volume={},
  number={},
  pages={2329-2333},
  doi={10.1109/PIERS.2017.8262141}}
```


## Installation
1. add the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) to your active MATLAB libraries.

2. run ```HELMHOLTZ_integral.m```

The code runs a sweep of experiments and stored results in your active MATLAB state. 