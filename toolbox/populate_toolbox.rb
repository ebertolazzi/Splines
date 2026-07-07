#!/usr/bin/env ruby

require "fileutils"

TOOLBOX_DIR = File.expand_path(__dir__)
ROOT_DIR    = File.expand_path("..", TOOLBOX_DIR)
DEPS_DIR    = File.expand_path("..", ROOT_DIR)

GC_DIR      = File.join(DEPS_DIR, "GC")
UTILS_DIR   = File.join(DEPS_DIR, "UtilsLite")
ROOTS_DIR   = File.join(DEPS_DIR, "quarticRootsFlocke")
JSON_DIR    = File.join(DEPS_DIR, "json")

SRC_DIR     = File.join(TOOLBOX_DIR, "src")
SRC_MEX_DIR = File.join(TOOLBOX_DIR, "src_mex")
BIN_DIR     = File.join(TOOLBOX_DIR, "bin")

def require_path(path)
  return path if File.exist?(path)

  raise "Missing required path: #{path}"
end

def copy_tree(src, dst)
  FileUtils.mkdir_p(dst)
  FileUtils.cp_r("#{src}/.", dst)
end

FileUtils.rm_rf(SRC_DIR)
Dir.glob(File.join(BIN_DIR, "*.mex*")).each { |file| File.delete(file) }

copy_tree(require_path(File.join(ROOT_DIR, "src")), SRC_DIR)
copy_tree(require_path(File.join(ROOT_DIR, "include")), SRC_DIR)
copy_tree(require_path(File.join(ROOTS_DIR, "src")), SRC_DIR)
copy_tree(require_path(File.join(UTILS_DIR, "src")), SRC_DIR)
copy_tree(require_path(File.join(GC_DIR, "src")), SRC_DIR)
copy_tree(require_path(File.join(GC_DIR, "include")), SRC_DIR)
copy_tree(require_path(File.join(JSON_DIR, "include")), SRC_DIR)

FileUtils.cp(
  require_path(File.join(GC_DIR, "matlab", "GenericContainerInterface_matlab.cc")),
  File.join(SRC_DIR, "GenericContainerInterface_matlab.cc")
)
FileUtils.cp(
  require_path(File.join(GC_DIR, "matlab", "GenericContainerInterface_matlab.cc")),
  File.join(SRC_MEX_DIR, "GenericContainerInterface_matlab.cc")
)

[
  File.join(SRC_DIR, "GenericContainerInterface_matlab.cc"),
  File.join(SRC_MEX_DIR, "GenericContainerInterface_matlab.cc"),
].each do |path|
  content = File.read(path)
  updated = content.gsub("GC_ASSERT(", "GC_assert(")
  File.write(path, updated) unless updated == content
end

# MATLAB toolbox version stays self-contained and does not vendor Eigen.
FileUtils.rm_rf(File.join(SRC_DIR, "Eigen"))
FileUtils.rm_rf(File.join(SRC_DIR, "unsupported"))

# Legacy cleanup kept from the original script.
[
  "Utils_Poly.cc",
  "Utils_GG2D.cc",
  "Utils_HJPatternSearch.cc",
  "Utils_NelderMead.cc",
  "Utils_nonlinear_system_tests.cc",
].each do |basename|
  FileUtils.rm_f(File.join(SRC_DIR, basename))
end

FileUtils.cp(require_path(File.join(ROOT_DIR, "license.txt")), File.join(TOOLBOX_DIR, "license.txt"))
