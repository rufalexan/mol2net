# mol2net
[![DOI](https://zenodo.org/badge/529230671.svg)](https://zenodo.org/badge/latestdoi/529230671)

Networks from a list of molecules

![fig_method](https://user-images.githubusercontent.com/112173397/186894549-131b817f-b398-404f-83e6-f362415c16d7.png)

## How to use this code

```python
python 1_2_mol2net.py  # Creates networks from a list of molecules
python 3_netAna.py  # Analyzes the network via networkX
```
Run the two scripts either separately, or automatically via the bash script:
```
mol2net.sh
```

## Credits
If you use this code for your work, please cite the corresponding paper
```
@article{ruf2022network,
author = {Ruf, Alexander and Danger, Grégoire},
title = {Network Analysis Reveals Spatial Clustering and Annotation of Complex Chemical Spaces: Application to Astrochemistry},
journal = {Analytical Chemistry},
volume = {94},
number = {41},
pages = {14135-14142},
year = {2022},
doi = {10.1021/acs.analchem.2c01271},
note ={PMID: 36209417},
URL = {https://doi.org/10.1021/acs.analchem.2c01271},
eprint = {https://doi.org/10.1021/acs.analchem.2c01271},
keywords={year2022,corresponding,first},
}
```
, and this github repository or the corresponding Zenodo DOI.
```
@software{rufalexan_2022_7025094,
  author       = {rufalexan},
  title        = {rufalexan/mol2net: v0.1.0},
  month        = aug,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.7025094},
  url          = {https://doi.org/10.5281/zenodo.7025094}
}
```

## Licence
This code is licensed under the MIT License.

