from __future__ import print_function
from _mbd_backend import ffi, lib as mbd
from mpi4py import MPI
import numpy as np
from itertools import chain
import json
import os


bohr = mbd.bohr
mbd.my_task = MPI.COMM_WORLD.Get_rank()
mbd.n_tasks = MPI.COMM_WORLD.Get_size()


def make_iter(obj, siz):
    for i in range(siz):
        yield obj[i]


def get_ndarray(pointer, shape, dtype=None):
    return np.ndarray(
        buffer=ffi.buffer(pointer, np.prod(shape)*8),
        shape=shape,
        dtype=dtype,
        order='F'
    )


class CStructWrapper(object):
    def __str__(self):
        s = ''
        for key, val in chain(
            ((attr, getattr(self, attr)) for attr in self._cattrs),
            ((attr, getattr(self, attr)) for attr in self._farrays),
            self.__dict__.items()
        ):
            if key.startswith('_'):
                continue
            strval = str(val)
            if '\n' in strval:
                s += '{0}:\n{1}\n'.format(key, strval)
            else:
                s += '{0}: {1}\n'.format(key, strval)
        return s.rstrip()

    def __getattr__(self, attr):
        try:
            return self.__dict__[attr]
        except KeyError:
            pass
        if attr in self._farrays:
            return get_ndarray(getattr(self._cobj, attr), **self._farrays[attr])
        if attr in self._cattrs:
            obj = getattr(self._cobj, attr)
            if isinstance(obj, ffi.CData) and ffi.typeof(obj) is ffi.typeof('char *'):
                return ffi.string(obj).decode()
            return obj
        raise AttributeError(
            "'{0}' object has no attribute '{1}'"
            .format(self.__class__.__name__, attr)
        )

    def __setattr__(self, attr, val):
        if attr in self.__dict__:
            raise AttributeError(
                "Cannot set attribute '{0}' of '{1}' object"
                .format(attr, self.__class__.__name__)
            )
        if attr in self._cattrs:
            return setattr(self._cobj, attr, val)
        self.__dict__[attr] = val


class Context(CStructWrapper):
    _cattrs = [
        'do_rpa', 'do_reciprocal', 'do_ewald', 'ts_energy_accuracy',
        'ts_shell_thickness', 'dipole_low_dim_cutoff', 'mayer_scaling',
        'ewald_real_cutoff_scaling', 'ewald_rec_cutoff_scaling'
    ]
    _farrays = []

    def __init__(self, n_grid=20):
        self._cobj = ffi.new('struct Context *')
        mbd._init_ctx(self._cobj, n_grid)
        self._farrays = {
            'freq_grid': {'shape': (n_grid,)},
            'freq_grid_w': {'shape': (n_grid,)}
        }

    def __del__(self):
        mbd._destroy_grid(self._cobj)


class Damping(CStructWrapper):
    _cattrs = [
        'label', 'd', 's_R', 'a', 'beta',
    ]

    def __init__(self, label, method, xc):
        self._cobj = ffi.new('struct Damping *')
        self._label = ffi.new('char[]', label.encode())
        self._cobj.label = self._label
        mbd._set_damping_parameters(
            self._cobj,
            ffi.new('char[]', method.encode()),
            ffi.new('char[]', xc.encode())
        )

    @property
    def _farrays(self):
        n = self._cobj.n_atoms
        return {
            'alpha_0': {'shape': (n,)},
            'C6': {'shape': (n,)},
            'R_vdw': {'shape': (n,)},
            'overlap': {'shape': (n, n)},
            'custom': {'shape': (n, n)},
            'custom_potential': {'shape': (n, n, n, n)}
        }


class Geometry(CStructWrapper):
    _cattrs = ['is_periodic']
    _farrays = {
        'supercell': {'shape': (3,), 'dtype': int},
        'lattice_vectors': {'shape': (3, 3)},
        'vacuum_axes': {'shape': (3,), 'dtype': bool}
    }

    with open(__file__[:-1] if __file__.endswith('.pyc') else __file__) as f:
        species_data = json.loads(
            next(l for l in f if l.startswith('# species_data:')).split(' ', 2)[-1]
        )
        free_atoms = json.loads(
            next(l for l in f if l.startswith('# free_atoms:')).split(' ', 2)[-1]
        )

    def __init__(self, species, coords, lattice=None):
        self._cobj = ffi.new('struct Geometry *')
        self.species = species
        self.coords = np.array(coords, order='F')
        self._cobj.coords = ffi.cast("double *", self.coords.ctypes.data)
        self._cobj.n_atoms = len(self)
        self.supercell[:] = 0
        self.lattice_vectors[:] = 0
        self.vacuum_axes[:] = False
        if lattice:
            self.lattice_vectors[:] = lattice
            self.is_periodic = True

    def __iter__(self):
        for specie, coord in zip(self.species, self.coords):
            yield specie, coord

    def __len__(self):
        return len(self.species)

    @classmethod
    def from_file(cls, f, fmt):
        if fmt == 'xyz':
            n = int(f.readline())
            f.readline()
            species = []
            coords = []
            for _ in range(n):
                l = f.readline().split()
                species.append(l[0])
                coords.append([float(x) for x in l[1:4]])
            return cls(species, coords)
        elif fmt == 'aims':
            species = []
            coords = []
            lattice = []
            while True:
                l = f.readline()
                if not l:
                    break
                l = l.strip()
                if not l or l.startswith('#'):
                    continue
                tokens = l.split()
                item = tokens[0]
                if item == 'atom':
                    species.append(tokens[4])
                    coords.append([float(x) for x in tokens[1:4]])
                elif item == 'lattice_vector':
                    lattice.append([float(x) for x in tokens[1:4]])
            if lattice:
                assert len(lattice) == 3
            return cls(species, coords, lattice)
        else:
            raise ValueError("Unknown format: '{0}'".format(fmt))

    @classmethod
    def from_path(cls, path, fmt=None):
        filename = os.path.basename(path)
        if not fmt:
            if filename == 'geometry.in':
                fmt = 'aims'
        if not fmt:
            _, ext = os.path.splitext(filename)
            fmt = {
                '.xyz': 'xyz',
                '.aims': 'aims'
            }.get(ext)
        if not fmt:
            raise ValueError("Unknown file format: '{0}'".format(filename))
        with open(path) as f:
            return cls.from_file(f, fmt)

    @classmethod
    def get_property(cls, idx, name):
        if isinstance(idx, int):
            select, value = 'number', str(idx)
        else:
            select, value = 'symbol', idx
        row = next(row for row in cls.species_data if row[select] == value)
        value = row[name]
        try:
            value = int(value)
        except ValueError:
            pass
        else:
            return value
        try:
            value = float(value)
        except ValueError:
            pass
        else:
            return value
        return value

    def get_free_atoms(self):
        return tuple(np.array([
            (sp['alpha_0'], sp['C6'], sp['R_vdw'])
            for sp in [self.free_atoms[sp] for sp in self.species]
        ]).T)


# species_data: [{"covalent_radius":"0.38","symbol":"H","number":"1","vdw_radius":"1.2","name":"hydrogen","ionization_energy":"13.5984","mass":"1.0079"},{"covalent_radius":"0.32","symbol":"He","number":"2","vdw_radius":"1.4","name":"helium","ionization_energy":"24.5874","mass":"4.0026"},{"covalent_radius":"1.34","symbol":"Li","number":"3","vdw_radius":"1.82","name":"lithium","ionization_energy":"5.3917","mass":"6.941"},{"covalent_radius":"0.9","symbol":"Be","number":"4","vdw_radius":"1.53","name":"beryllium","ionization_energy":"9.3227","mass":"9.0122"},{"covalent_radius":"0.82","symbol":"B","number":"5","vdw_radius":"1.92","name":"boron","ionization_energy":"8.298","mass":"10.811"},{"covalent_radius":"0.77","symbol":"C","number":"6","vdw_radius":"1.7","name":"carbon","ionization_energy":"11.2603","mass":"12.0107"},{"covalent_radius":"0.75","symbol":"N","number":"7","vdw_radius":"1.55","name":"nitrogen","ionization_energy":"14.5341","mass":"14.0067"},{"covalent_radius":"0.73","symbol":"O","number":"8","vdw_radius":"1.52","name":"oxygen","ionization_energy":"13.6181","mass":"15.9994"},{"covalent_radius":"0.71","symbol":"F","number":"9","vdw_radius":"1.47","name":"fluorine","ionization_energy":"17.4228","mass":"18.9984"},{"covalent_radius":"0.69","symbol":"Ne","number":"10","vdw_radius":"1.54","name":"neon","ionization_energy":"21.5645","mass":"20.1797"},{"covalent_radius":"1.54","symbol":"Na","number":"11","vdw_radius":"2.27","name":"sodium","ionization_energy":"5.1391","mass":"22.9897"},{"covalent_radius":"1.3","symbol":"Mg","number":"12","vdw_radius":"1.73","name":"magnesium","ionization_energy":"7.6462","mass":"24.305"},{"covalent_radius":"1.18","symbol":"Al","number":"13","vdw_radius":"1.84","name":"aluminium","ionization_energy":"5.9858","mass":"26.9815"},{"covalent_radius":"1.11","symbol":"Si","number":"14","vdw_radius":"2.1","name":"silicon","ionization_energy":"8.1517","mass":"28.0855"},{"covalent_radius":"1.06","symbol":"P","number":"15","vdw_radius":"1.8","name":"phosphorus","ionization_energy":"10.4867","mass":"30.9738"},{"covalent_radius":"1.02","symbol":"S","number":"16","vdw_radius":"1.8","name":"sulfur","ionization_energy":"10.36","mass":"32.065"},{"covalent_radius":"0.99","symbol":"Cl","number":"17","vdw_radius":"1.75","name":"chlorine","ionization_energy":"12.9676","mass":"35.453"},{"covalent_radius":"0.97","symbol":"Ar","number":"18","vdw_radius":"1.88","name":"argon","ionization_energy":"15.7596","mass":"39.948"},{"covalent_radius":"1.96","symbol":"K","number":"19","vdw_radius":"2.75","name":"potassium","ionization_energy":"4.3407","mass":"39.0983"},{"covalent_radius":"1.74","symbol":"Ca","number":"20","vdw_radius":"2.31","name":"calcium","ionization_energy":"6.1132","mass":"40.078"},{"covalent_radius":"1.44","symbol":"Sc","number":"21","vdw_radius":"2.11","name":"scandium","ionization_energy":"6.5615","mass":"44.9559"},{"covalent_radius":"1.36","symbol":"Ti","number":"22","vdw_radius":"","name":"titanium","ionization_energy":"6.8281","mass":"47.867"},{"covalent_radius":"1.25","symbol":"V","number":"23","vdw_radius":"","name":"vanadium","ionization_energy":"6.7462","mass":"50.9415"},{"covalent_radius":"1.27","symbol":"Cr","number":"24","vdw_radius":"","name":"chromium","ionization_energy":"6.7665","mass":"51.9961"},{"covalent_radius":"1.39","symbol":"Mn","number":"25","vdw_radius":"","name":"manganese","ionization_energy":"7.434","mass":"54.938"},{"covalent_radius":"1.25","symbol":"Fe","number":"26","vdw_radius":"","name":"iron","ionization_energy":"7.9024","mass":"55.845"},{"covalent_radius":"1.26","symbol":"Co","number":"27","vdw_radius":"","name":"cobalt","ionization_energy":"7.881","mass":"58.9332"},{"covalent_radius":"1.21","symbol":"Ni","number":"28","vdw_radius":"1.63","name":"nickel","ionization_energy":"7.6398","mass":"58.6934"},{"covalent_radius":"1.38","symbol":"Cu","number":"29","vdw_radius":"1.4","name":"copper","ionization_energy":"7.7264","mass":"63.546"},{"covalent_radius":"1.31","symbol":"Zn","number":"30","vdw_radius":"1.39","name":"zinc","ionization_energy":"9.3942","mass":"65.39"},{"covalent_radius":"1.26","symbol":"Ga","number":"31","vdw_radius":"1.87","name":"gallium","ionization_energy":"5.9993","mass":"69.723"},{"covalent_radius":"1.22","symbol":"Ge","number":"32","vdw_radius":"2.11","name":"germanium","ionization_energy":"7.8994","mass":"72.64"},{"covalent_radius":"1.19","symbol":"As","number":"33","vdw_radius":"1.85","name":"arsenic","ionization_energy":"9.7886","mass":"74.9216"},{"covalent_radius":"1.16","symbol":"Se","number":"34","vdw_radius":"1.9","name":"selenium","ionization_energy":"9.7524","mass":"78.96"},{"covalent_radius":"1.14","symbol":"Br","number":"35","vdw_radius":"1.85","name":"bromine","ionization_energy":"11.8138","mass":"79.904"},{"covalent_radius":"1.1","symbol":"Kr","number":"36","vdw_radius":"2.02","name":"krypton","ionization_energy":"13.9996","mass":"83.8"},{"covalent_radius":"2.11","symbol":"Rb","number":"37","vdw_radius":"3.03","name":"rubidium","ionization_energy":"4.1771","mass":"85.4678"},{"covalent_radius":"1.92","symbol":"Sr","number":"38","vdw_radius":"2.49","name":"strontium","ionization_energy":"5.6949","mass":"87.62"},{"covalent_radius":"1.62","symbol":"Y","number":"39","vdw_radius":"","name":"yttrium","ionization_energy":"6.2173","mass":"88.9059"},{"covalent_radius":"1.48","symbol":"Zr","number":"40","vdw_radius":"","name":"zirconium","ionization_energy":"6.6339","mass":"91.224"},{"covalent_radius":"1.37","symbol":"Nb","number":"41","vdw_radius":"","name":"niobium","ionization_energy":"6.7589","mass":"92.9064"},{"covalent_radius":"1.45","symbol":"Mo","number":"42","vdw_radius":"","name":"molybdenum","ionization_energy":"7.0924","mass":"95.94"},{"covalent_radius":"1.56","symbol":"Tc","number":"43","vdw_radius":"","name":"technetium","ionization_energy":"7.28","mass":"98"},{"covalent_radius":"1.26","symbol":"Ru","number":"44","vdw_radius":"","name":"ruthenium","ionization_energy":"7.3605","mass":"101.07"},{"covalent_radius":"1.35","symbol":"Rh","number":"45","vdw_radius":"","name":"rhodium","ionization_energy":"7.4589","mass":"102.9055"},{"covalent_radius":"1.31","symbol":"Pd","number":"46","vdw_radius":"1.63","name":"palladium","ionization_energy":"8.3369","mass":"106.42"},{"covalent_radius":"1.53","symbol":"Ag","number":"47","vdw_radius":"1.72","name":"silver","ionization_energy":"7.5762","mass":"107.8682"},{"covalent_radius":"1.48","symbol":"Cd","number":"48","vdw_radius":"1.58","name":"cadmium","ionization_energy":"8.9938","mass":"112.411"},{"covalent_radius":"1.44","symbol":"In","number":"49","vdw_radius":"1.93","name":"indium","ionization_energy":"5.7864","mass":"114.818"},{"covalent_radius":"1.41","symbol":"Sn","number":"50","vdw_radius":"2.17","name":"tin","ionization_energy":"7.3439","mass":"118.71"},{"covalent_radius":"1.38","symbol":"Sb","number":"51","vdw_radius":"2.06","name":"antimony","ionization_energy":"8.6084","mass":"121.76"},{"covalent_radius":"1.35","symbol":"Te","number":"52","vdw_radius":"2.06","name":"tellurium","ionization_energy":"9.0096","mass":"127.6"},{"covalent_radius":"1.33","symbol":"I","number":"53","vdw_radius":"1.98","name":"iodine","ionization_energy":"10.4513","mass":"126.9045"},{"covalent_radius":"1.3","symbol":"Xe","number":"54","vdw_radius":"2.16","name":"xenon","ionization_energy":"12.1298","mass":"131.293"},{"covalent_radius":"2.25","symbol":"Cs","number":"55","vdw_radius":"3.43","name":"caesium","ionization_energy":"3.8939","mass":"132.9055"},{"covalent_radius":"1.98","symbol":"Ba","number":"56","vdw_radius":"2.68","name":"barium","ionization_energy":"5.2117","mass":"137.327"},{"covalent_radius":"1.69","symbol":"La","number":"57","vdw_radius":"","name":"lanthanum","ionization_energy":"5.5769","mass":"138.9055"},{"covalent_radius":"","symbol":"Ce","number":"58","vdw_radius":"","name":"cerium","ionization_energy":"5.5387","mass":"140.116"},{"covalent_radius":"","symbol":"Pr","number":"59","vdw_radius":"","name":"praseodymium","ionization_energy":"5.473","mass":"140.9077"},{"covalent_radius":"","symbol":"Nd","number":"60","vdw_radius":"","name":"neodymium","ionization_energy":"5.525","mass":"144.24"},{"covalent_radius":"","symbol":"Pm","number":"61","vdw_radius":"","name":"promethium","ionization_energy":"5.582","mass":"145"},{"covalent_radius":"","symbol":"Sm","number":"62","vdw_radius":"","name":"samarium","ionization_energy":"5.6437","mass":"150.36"},{"covalent_radius":"","symbol":"Eu","number":"63","vdw_radius":"","name":"europium","ionization_energy":"5.6704","mass":"151.964"},{"covalent_radius":"","symbol":"Gd","number":"64","vdw_radius":"","name":"gadolinium","ionization_energy":"6.1501","mass":"157.25"},{"covalent_radius":"","symbol":"Tb","number":"65","vdw_radius":"","name":"terbium","ionization_energy":"5.8638","mass":"158.9253"},{"covalent_radius":"","symbol":"Dy","number":"66","vdw_radius":"","name":"dysprosium","ionization_energy":"5.9389","mass":"162.5"},{"covalent_radius":"","symbol":"Ho","number":"67","vdw_radius":"","name":"holmium","ionization_energy":"6.0215","mass":"164.9303"},{"covalent_radius":"","symbol":"Er","number":"68","vdw_radius":"","name":"erbium","ionization_energy":"6.1077","mass":"167.259"},{"covalent_radius":"","symbol":"Tm","number":"69","vdw_radius":"","name":"thulium","ionization_energy":"6.1843","mass":"168.9342"},{"covalent_radius":"","symbol":"Yb","number":"70","vdw_radius":"","name":"ytterbium","ionization_energy":"6.2542","mass":"173.04"},{"covalent_radius":"1.6","symbol":"Lu","number":"71","vdw_radius":"","name":"lutetium","ionization_energy":"5.4259","mass":"174.967"},{"covalent_radius":"1.5","symbol":"Hf","number":"72","vdw_radius":"","name":"hafnium","ionization_energy":"6.8251","mass":"178.49"},{"covalent_radius":"1.38","symbol":"Ta","number":"73","vdw_radius":"","name":"tantalum","ionization_energy":"7.5496","mass":"180.9479"},{"covalent_radius":"1.46","symbol":"W","number":"74","vdw_radius":"","name":"tungsten","ionization_energy":"7.864","mass":"183.84"},{"covalent_radius":"1.59","symbol":"Re","number":"75","vdw_radius":"","name":"rhenium","ionization_energy":"7.8335","mass":"186.207"},{"covalent_radius":"1.28","symbol":"Os","number":"76","vdw_radius":"","name":"osmium","ionization_energy":"8.4382","mass":"190.23"},{"covalent_radius":"1.37","symbol":"Ir","number":"77","vdw_radius":"","name":"iridium","ionization_energy":"8.967","mass":"192.217"},{"covalent_radius":"1.28","symbol":"Pt","number":"78","vdw_radius":"1.75","name":"platinum","ionization_energy":"8.9587","mass":"195.078"},{"covalent_radius":"1.44","symbol":"Au","number":"79","vdw_radius":"1.66","name":"gold","ionization_energy":"9.2255","mass":"196.9665"},{"covalent_radius":"1.49","symbol":"Hg","number":"80","vdw_radius":"1.55","name":"mercury","ionization_energy":"10.4375","mass":"200.59"},{"covalent_radius":"1.48","symbol":"Tl","number":"81","vdw_radius":"1.96","name":"thallium","ionization_energy":"6.1082","mass":"204.3833"},{"covalent_radius":"1.47","symbol":"Pb","number":"82","vdw_radius":"2.02","name":"lead","ionization_energy":"7.4167","mass":"207.2"},{"covalent_radius":"1.46","symbol":"Bi","number":"83","vdw_radius":"2.07","name":"bismuth","ionization_energy":"7.2856","mass":"208.9804"},{"covalent_radius":"","symbol":"Po","number":"84","vdw_radius":"1.97","name":"polonium","ionization_energy":"8.417","mass":"209"},{"covalent_radius":"","symbol":"At","number":"85","vdw_radius":"2.02","name":"astatine","ionization_energy":"9.3","mass":"210"},{"covalent_radius":"1.45","symbol":"Rn","number":"86","vdw_radius":"2.2","name":"radon","ionization_energy":"10.7485","mass":"222"},{"covalent_radius":"","symbol":"Fr","number":"87","vdw_radius":"3.48","name":"francium","ionization_energy":"4.0727","mass":"223"},{"covalent_radius":"","symbol":"Ra","number":"88","vdw_radius":"2.83","name":"radium","ionization_energy":"5.2784","mass":"226"},{"covalent_radius":"","symbol":"Ac","number":"89","vdw_radius":"","name":"actinium","ionization_energy":"5.17","mass":"227"},{"covalent_radius":"","symbol":"Th","number":"90","vdw_radius":"","name":"thorium","ionization_energy":"6.3067","mass":"232.0381"},{"covalent_radius":"","symbol":"Pa","number":"91","vdw_radius":"","name":"protactinium","ionization_energy":"5.89","mass":"231.0359"},{"covalent_radius":"","symbol":"U","number":"92","vdw_radius":"1.86","name":"uranium","ionization_energy":"6.1941","mass":"238.0289"}]
# free_atoms: {"H":{"alpha_0":4.500000,"C6":6.500000,"R_vdw":3.100000},"He":{"alpha_0":1.380000,"C6":1.460000,"R_vdw":2.650000},"Li":{"alpha_0":164.200000,"C6":1387.000000,"R_vdw":4.160000},"Be":{"alpha_0":38.000000,"C6":214.000000,"R_vdw":4.170000},"B":{"alpha_0":21.000000,"C6":99.500000,"R_vdw":3.890000},"C":{"alpha_0":12.000000,"C6":46.600000,"R_vdw":3.590000},"N":{"alpha_0":7.400000,"C6":24.200000,"R_vdw":3.340000},"O":{"alpha_0":5.400000,"C6":15.600000,"R_vdw":3.190000},"F":{"alpha_0":3.800000,"C6":9.520000,"R_vdw":3.040000},"Ne":{"alpha_0":2.670000,"C6":6.380000,"R_vdw":2.910000},"Na":{"alpha_0":162.700000,"C6":1556.000000,"R_vdw":3.730000},"Mg":{"alpha_0":71.000000,"C6":627.000000,"R_vdw":4.270000},"Al":{"alpha_0":60.000000,"C6":528.000000,"R_vdw":4.330000},"Si":{"alpha_0":37.000000,"C6":305.000000,"R_vdw":4.200000},"P":{"alpha_0":25.000000,"C6":185.000000,"R_vdw":4.010000},"S":{"alpha_0":19.600000,"C6":134.000000,"R_vdw":3.860000},"Cl":{"alpha_0":15.000000,"C6":94.600000,"R_vdw":3.710000},"Ar":{"alpha_0":11.100000,"C6":64.300000,"R_vdw":3.550000},"K":{"alpha_0":292.900000,"C6":3897.000000,"R_vdw":3.710000},"Ca":{"alpha_0":160.000000,"C6":2221.000000,"R_vdw":4.650000},"Sc":{"alpha_0":120.000000,"C6":1383.000000,"R_vdw":4.590000},"Ti":{"alpha_0":98.000000,"C6":1044.000000,"R_vdw":4.510000},"V":{"alpha_0":84.000000,"C6":832.000000,"R_vdw":4.440000},"Cr":{"alpha_0":78.000000,"C6":602.000000,"R_vdw":3.990000},"Mn":{"alpha_0":63.000000,"C6":552.000000,"R_vdw":3.970000},"Fe":{"alpha_0":56.000000,"C6":482.000000,"R_vdw":4.230000},"Co":{"alpha_0":50.000000,"C6":408.000000,"R_vdw":4.180000},"Ni":{"alpha_0":48.000000,"C6":373.000000,"R_vdw":3.820000},"Cu":{"alpha_0":42.000000,"C6":253.000000,"R_vdw":3.760000},"Zn":{"alpha_0":40.000000,"C6":284.000000,"R_vdw":4.020000},"Ga":{"alpha_0":60.000000,"C6":498.000000,"R_vdw":4.190000},"Ge":{"alpha_0":41.000000,"C6":354.000000,"R_vdw":4.200000},"As":{"alpha_0":29.000000,"C6":246.000000,"R_vdw":4.110000},"Se":{"alpha_0":25.000000,"C6":210.000000,"R_vdw":4.040000},"Br":{"alpha_0":20.000000,"C6":162.000000,"R_vdw":3.930000},"Kr":{"alpha_0":16.800000,"C6":129.600000,"R_vdw":3.820000},"Rb":{"alpha_0":319.200000,"C6":4691.000000,"R_vdw":3.720000},"Sr":{"alpha_0":199.000000,"C6":3170.000000,"R_vdw":4.540000},"Y":{"alpha_0":126.7370,"C6":1968.580,"R_vdw":4.81510},"Zr":{"alpha_0":119.97,"C6":1677.91,"R_vdw":4.53},"Nb":{"alpha_0":101.603,"C6":1263.61,"R_vdw":4.2365},"Mo":{"alpha_0":88.4225785,"C6":1028.73,"R_vdw":4.099},"Tc":{"alpha_0":80.083,"C6":1390.87,"R_vdw":4.076},"Ru":{"alpha_0":65.8950,"C6":609.754,"R_vdw":3.99530},"Rh":{"alpha_0":56.1,"C6":469.0,"R_vdw":3.95},"Pd":{"alpha_0":23.680000,"C6":157.500000,"R_vdw":3.66000},"Ag":{"alpha_0":50.600000,"C6":339.000000,"R_vdw":3.820000},"Cd":{"alpha_0":39.7,"C6":452.0,"R_vdw":3.99},"In":{"alpha_0":70.22000,"C6":707.046000,"R_vdw":4.23198000},"Sn":{"alpha_0":55.9500,"C6":587.41700,"R_vdw":4.303000},"Sb":{"alpha_0":43.671970,"C6":459.322,"R_vdw":4.2760},"Te":{"alpha_0":37.65,"C6":396.0,"R_vdw":4.22},"I":{"alpha_0":35.000000,"C6":385.000000,"R_vdw":4.170000},"Xe":{"alpha_0":27.300000,"C6":285.900000,"R_vdw":4.080000},"Cs":{"alpha_0":427.12,"C6":6582.08,"R_vdw":3.78},"Ba":{"alpha_0":275.0,"C6":5727.0,"R_vdw":4.77},"Hf":{"alpha_0":99.52,"C6":1274.8,"R_vdw":4.21},"Ta":{"alpha_0":82.53,"C6":1019.92,"R_vdw":4.15},"W":{"alpha_0":71.041,"C6":847.93,"R_vdw":4.08},"Re":{"alpha_0":63.04,"C6":710.2,"R_vdw":4.02},"Os":{"alpha_0":55.055,"C6":596.67,"R_vdw":3.84},"Ir":{"alpha_0":42.51,"C6":359.1,"R_vdw":4.00},"Pt":{"alpha_0":39.68,"C6":347.1,"R_vdw":3.92},"Au":{"alpha_0":36.5,"C6":298.0,"R_vdw":3.86},"Hg":{"alpha_0":33.9,"C6":392.0,"R_vdw":3.98},"Tl":{"alpha_0":69.92,"C6":717.44,"R_vdw":3.91},"Pb":{"alpha_0":61.8,"C6":697.0,"R_vdw":4.31},"Bi":{"alpha_0":49.02,"C6":571.0,"R_vdw":4.32},"Po":{"alpha_0":45.013,"C6":530.92,"R_vdw":4.097},"At":{"alpha_0":38.93,"C6":457.53,"R_vdw":4.07},"Rn":{"alpha_0":33.54,"C6":390.63,"R_vdw":4.23}}
