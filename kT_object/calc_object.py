from ase.calculators.vasp import Vasp

class default_calc:
    def __init__(self):
        pass
    def load_default_vasp(self):
        self.vasp_params = {"kpts": (1,1,1),
                       "potim":0.5,
                       "encut":500,
                       "ispin":2,
                       "nsw":50,
                       "prec":'Normal',
                       "istart":1,
                       "isif":2,
                       "ismear":0,
                       "sigma":0.05,
                       "nelmin":4,
                       "nelmdl":-4,
                       "nwrite":1,
                       "icharg":2,
                       "lasph":True,
                       "ediff":1E-6,
                       "ediffg":-0.05, 
                       "ibrion":2,
                       "lcharg":False,
                       "lwave":False,
                       "laechg":False,
                       "voskown":1,
                       "algo":'Fast',
                       "lplane":True,
                       "lreal":'Auto',
                       "isym":0,
                       "xc":'PBE',
                       "lorbit":11,
                       "nupdown":-1,
                       "npar":4,
                       "nsim":4,
                       "ivdw":12
                       }
        return self.vasp_params
