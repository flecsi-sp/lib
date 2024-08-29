from spack.package import *

class FlecsiSP(CMakePackage, CudaPackage, ROCmPackage):
    """The FleCSI-SP library provides utilities for creating FleCSI specializations
    """

    homepage = "http://flecsi-sp.org/"
    git = "https://github.com/flecsi-sp/lib.git"
    maintainers("bergen")

    version("develop", branch="develop")
