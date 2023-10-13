# Towards Robust Hyperspectral Unmixing: Mixed Noise Modeling and Image-Domain Regularization

This is a demo code of the method proposed in the following reference:

K. Naganuma and S. Ono
``Towards Robust Hyperspectral Unmixing: Mixed Noise Modeling and Image-Domain Regularization''

Update history:
Octber 13, 2023: v1.0 

For more information, see the following 
- Project website: https://www.mdi.c.titech.ac.jp/publications/rhuidr
- Preprint paper: https://arxiv.org/abs/2302.08247

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
      title={Towards Robust Hyperspectral Unmixing: Mixed Noise Modeling and Image-Domain Regularization}, 
      author={Kazuki Naganuma and Shunsuke Ono},
      year={2023},
      eprint={2302.08247},
      archivePrefix={arXiv},
      primaryClass={eess.IV}
}
```
