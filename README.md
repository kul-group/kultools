# KulTools
Kultools is an interface to run atomistic simulations using [ASE](https://wiki.fysik.dtu.dk/ase/index.html) and [VASP](https://www.vasp.at/wiki/index.php/The_VASP_Manual). 

⚠️ Documentation under construction 



[![PyPI](https://img.shields.io/pypi/v/kultools)](https://pypi.org/project/kultools/) 
[![Python-Version](https://img.shields.io/badge/Python-3.6+-green)](https://github.com/kul-group/kultools)
[![Downloads](https://static.pepy.tech/badge/kultools/month)](https://pepy.tech/project/kultools)
[![GitHub contributors](https://img.shields.io/github/contributors/kul-group/kultools)](https://github.com/kul-group/kultools/graphs/contributors)


## Install

The package can be installed using [pip](https://pypi.org/project/kultools/). You would require python 3.6 or above

```bash
pip install kultools
```

## Usage

**Note:** Examples coming soon!

### Running a zeolite calculation using KulTools

```python
>>> from kul_tools import KulTools as KT
>>> from ase.build import molecule

>>> kt = KT(gamma_only=False,structure_type='zeo',is_stop_eligible=True)
KT: HPC= local
KT: VASP_GAMMA= False
KT: VASP_PP_PATH= local_vasp_pp
KT: VASP_COMMAND= local_vasp_std
>>> kt.set_calculation_type('opt')
>>> atoms = molecule('H2O')
>>> atoms.set_cell(8*np.identity(3))
>>> atoms.center()
>>> atoms.pbc=True
>>> 
>>> kt.set_structure(atoms)
>>> kt.set_overall_vasp_params({'gga':'RP'})
>>> atoms = kt.run()

```

## Task List

- [ ] Update code documentation
- [ ] Add examples
- [ ] Update Readme
- [ ] Setup GitHub Actions

## Support and contribution

If you like the project, please ⭑ the project! Also, feel free to report unexpected behaviour at [issue](https://github.com/kul-group/kultools/issues). 

Contributions are always welcome!! Please open a [PR](https://github.com/kul-group/kultools/pulls) with corresponding changes to documentation and tests.


## Side notes

- load vasp module
