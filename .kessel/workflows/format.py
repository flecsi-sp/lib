from pathlib import Path

from kessel.workflows import collapsed, environment
from kessel.workflows.base.cmake import CMake
from kessel.workflows.base.spack import BuildEnvironment


class Format(BuildEnvironment, CMake):
    steps = ["env", "configure", "check_format"]
    allow_lockfile_changes = True

    build_dir = environment(Path.cwd() / "build_format")
    spack_env = environment("flecsi-sp-format")
    project_spec = environment("flecsi-sp+format")

    def ci_message(self, args):
        return super().ci_message(args, post_alloc_init="source .gitlab/kessel.sh")

    @collapsed
    def configure(self, args):
        """Configure"""
        cmake_args = [
            self.define("ENABLE_WARNINGS", False),
            self.define("ENABLE_FORMAT", True),
        ]
        super().configure(args, cmake_args)

    def check_format(self, args):
        """Clang-Format"""
        super().build(args, targets=["format"])
