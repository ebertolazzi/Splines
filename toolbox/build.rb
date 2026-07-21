require "fileutils"

TOOLBOX_DIR = File.expand_path(__dir__)
ROOT_DIR    = File.expand_path("..", TOOLBOX_DIR)
BUILD_DIR   = File.join(TOOLBOX_DIR, "build")
DEPS_BUILD  = File.join(BUILD_DIR, "dependencies")
MEX_BUILD   = File.join(BUILD_DIR, "matlab")

def run!(*command)
  return if system(*command)

  abort "Command failed: #{command.join(' ')}"
end

FileUtils.rm_rf(BUILD_DIR)

# Configuring the main project resolves local checkouts or FetchContent
# dependencies and populates toolbox/src entirely through CMake.
run!(
  "cmake", "-S", ROOT_DIR, "-B", DEPS_BUILD,
  "-DSPLINES_POPULATE_TOOLBOX=ON",
  "-DSPLINES_UPDATE_3RDPARTY=OFF",
  "-DSPLINES_INSTALL=OFF",
  "-DSPLINES_BUILD_BENCHMARKS=OFF",
  "-DBUILD_TESTING=OFF"
)

run!("cmake", "-S", TOOLBOX_DIR, "-B", MEX_BUILD)
run!("cmake", "--build", MEX_BUILD, "--parallel")
