from spack.package import *


class FlecsiSp(CMakePackage):
    """The FleCSI-SP library provides utilities for creating FleCSI specializations
    """

    ############################################################################
    # Info
    ############################################################################

    homepage = "http://flecsi-sp.org/"
    git = "ssh://git@re-git.lanl.gov:10022/flecsi-sp/lib.git"
    maintainers("bergen")

    ############################################################################
    # Versions
    ############################################################################

    version("develop", branch="develop")

    ############################################################################
    # Variants
    ############################################################################

    variant("exodusii", default=True,
            description="Build with support for the ExodusII file format"
            )
    variant("x3d", default=True,
            description="Build with support for the X3D file format"
            )
    variant("format", default=False,
            description="add dependencies for format check")
    variant("documentation", default=False, description="Enable documentation")
    variant("zoltan", default=False, description="Enable domain decomposition via Zoltan")

    ############################################################################
    # Dependencies
    ############################################################################

    depends_on("flecsi@2.4.0:")

    depends_on("exodusii", when="+exodusii")

    depends_on("cmake@3.27:")
    depends_on("py-sphinx", when="+documentation", type="build")
    depends_on("py-sphinx-rtd-theme", when="+documentation", type="build")
    depends_on("doxygen", when="+documentation", type="build")
    depends_on("graphviz", when="+documentation", type="build")
    depends_on("llvm@18", when="+format", type="build")
    depends_on("parmetis@4.0.3:")
    depends_on("zoltan+parmetis+mpi", when="+zoltan")

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    ############################################################################
    # Build
    ############################################################################

    def cmake_args(self):
        spec = self.spec

        options = [
            self.define_from_variant("ENABLE_EXODUSII", "exodusii"),
            self.define_from_variant("ENABLE_X3D", "x3d"),
            self.define_from_variant("ENABLE_DOCUMENTATION", "documentation"),
            self.define_from_variant("ENABLE_ZOLTAN", "zoltan"),
            self.define("ENABLE_UNIT_TESTS", self.run_tests)
        ]

        return options
