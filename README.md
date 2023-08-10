# Robust_Constrained_Hyperspectral_Unmixing_Using_Reconstructed-Image_Regularization

This is the demo code of the method proposed in the following reference:

K. Naganuma, Y. Nagamatsu, and S. Ono
``Robust Constrained Hyperspectral Unmixing Using Reconstructed-Image Regularization.''

Update history:
Augast 7, 2023: v1.0 

For more information, see the following 
- the project website: https://www.mdi.c.titech.ac.jp/publications/rchu
- the preprint paper: https://arxiv.org/abs/2302.08247

# How to use
Run `main.m`

# Contents
#### Data
We use the HYperspectral Data Retrieval and Analysis (HYDRA) toolbox to generate a synthetic HS image.
HYDRA can be obtained at https://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Imagery_Synthesis_tools_for_MATLAB.

#### Spectral library
To make an endmember library, we use spectral signatures from the U.S. Geological Survey (USGS) Spectral Library accessed at https://www.usgs.gov/programs/usgs-library.

# Reference
If you use this code, please cite the following paper:

```
@misc{naganuma2023robust,
      title={Robust Constrained Hyperspectral Unmixing Using Reconstructed-Image Regularization}, 
      author={Kazuki Naganuma and Yuki Nagamatsu and Shunsuke Ono},
      year={2023},
      eprint={2302.08247},
      archivePrefix={arXiv},
      primaryClass={eess.IV}
}
```
