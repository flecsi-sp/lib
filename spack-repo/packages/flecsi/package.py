from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI.
    """
    version("2.4-devel", commit="0c3949e136fe2c9739829e9f606ccaab7cc3ac6e")
