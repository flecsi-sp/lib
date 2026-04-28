from kessel.workflows import environment
from kessel.workflows.base.cmake import CMake
from kessel.workflows.base.spack import BuildEnvironment


class Default(BuildEnvironment, CMake):
    steps = ["env", "configure", "build", "test", "install"]
    allow_lockfile_changes = True

    project_spec = environment("flecsi-sp")

    def ci_message(self, args):
        return super().ci_message(args, post_alloc_init="source .gitlab/kessel.sh")

    def build(self, args):
        """Build"""
        cmake_args = [
            self.define("CMAKE_BUILD_TYPE", "Debug"),
            self.define("ENABLE_UNIT_TESTS", True),
        ]
        super().build(args, cmake_args)
