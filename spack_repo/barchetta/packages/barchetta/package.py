from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *

class Barchetta(CMakePackage):
    """
    Barchetta package
    """

    ############################################################################
    # Info
    ############################################################################

    homepage = "https://re-git.lanl.gov/rush/barchetta"
    git = "ssh://git@asc-git.lanl.gov:10022/rush/barchetta.git"
    maintainers("ayaguelopez")

    ############################################################################
    # Versions
    ############################################################################

    version("develop", branch="develop", submodules=True)

    ############################################################################
    # Variants
    ############################################################################

    variant("precision", default="double", values=("float", "double"),
        description="Select the precision", multi=False)
    variant("spiner", default=False, description="enable support tabular EOS")
    variant("format", default=False,
            description="add dependencies for format check")
    variant("documentation", default=False,
            description="Build radio documentation")

    ############################################################################
    # Dependencies
    ############################################################################

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on("flecsi-sp")
    depends_on("flecsi@2.4.0:")
    depends_on("py-pybind11")
    depends_on("mfem@:4.8+shared")
    depends_on("singularity-eos@main", when="~spiner")
    depends_on("singularity-eos@main+spiner+eospac+hdf5 build_extra=sesame",
               when="+spiner")
    depends_on("py-sphinx", when="+documentation", type="build")
    depends_on("py-sphinx-rtd-theme", when="+documentation", type="build")
    depends_on("py-breathe", when="+documentation", type="build")
    depends_on("py-sphinxcontrib-tikz", when="+documentation", type="build")
    depends_on("llvm@18", when="+format", type="build")

    ############################################################################
    # Build
    ############################################################################

    def cmake_args(self):

        options = [
            self.define_from_variant("CONFIG_PRECISION", "precision"),
            self.define("ENABLE_UNIT_TESTS", self.run_tests),
            self.define_from_variant("ENABLE_DOCUMENTATION", "documentation")
        ]

        return options
