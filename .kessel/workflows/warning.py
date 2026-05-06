from pathlib import Path

from kessel.workflows import collapsed, environment
from kessel.workflows.base.cmake import CMake
from kessel.workflows.base.spack import BuildEnvironment


class Warning(BuildEnvironment, CMake):
    steps = ["env", "configure", "check_warnings"]
    allow_lockfile_changes = True

    build_dir = environment(Path.cwd() / "build_warnings")
    spack_env = environment("flecsi-sp-warning")
    project_spec = environment("flecsi-sp")

    def ci_message(self, args):
        return super().ci_message(args, post_alloc_init="source .gitlab/kessel.sh")

    @collapsed
    def configure(self, args):
        """Configure"""
        cmake_args = [
            self.define("ENABLE_WARNINGS", True),
            self.define("WARNINGS_FATAL", True),
        ]
        super().configure(args, cmake_args)

    def check_warnings(self, args):
        """Check-Warnings"""
        super().build(args)
