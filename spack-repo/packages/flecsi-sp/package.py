from spack.package import *

class FlecsiSp(CMakePackage):
    """The FleCSI-SP library provides utilities for creating FleCSI specializations
    """

    homepage = "http://flecsi-sp.org/"
    git = "ssh://git@re-git.lanl.gov:10022/flecsi-sp/lib.git"
    maintainers("bergen")

    version("develop", branch="develop")

    variant("exodusii", default=True,
            description="Build with support for the ExodusII file format"
    )
    variant("x3d", default=True,
            description="Build with support for the X3D file format"
    )
    variant("documentation", default=False, description="Enable documentation")

    depends_on("flecsi@2.3:")

    depends_on("exodusii", when="+exodusii")

    depends_on("py-sphinx", when="+documentation", type="build")
    depends_on("py-sphinx-rtd-theme", when="+documentation", type="build")
    depends_on("doxygen", when="+documentation", type="build")
    depends_on("graphviz", when="+documentation", type="build")

    def cmake_args(self):
        spec = self.spec

        options = [
            self.define_from_variant("ENABLE_EXODUSII", "exodusii"),
            self.define_from_variant("ENABLE_X3D", "x3d"),
            self.define_from_variant("ENABLE_DOCUMENTATION", "doc"),
            self.define("ENABLE_UNIT_TESTS", self.run_tests)
        ]

        return options
