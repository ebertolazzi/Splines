#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"
lst = Dir["../src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../GC/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../GC/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

FileUtils.cp "../GC/src_matlab_interface/GenericContainerMatlabInterface.cc",
             "./src/GenericContainerMatlabInterface.cc"
FileUtils.cp "../GC/src_matlab_interface/GenericContainerMatlabInterface.hh", 
             "./src/GenericContainerMatlabInterface.hh"

FileUtils.rm_rf   "lib"
FileUtils.mkdir_p "lib"
lst = Dir["../lib_matlab/*.m"]
lst.each do |filename|
  FileUtils.cp filename, "./lib/" + File.basename(filename);
end

FileUtils.rm_rf   "src_mex"
FileUtils.mkdir_p "src_mex"
lst = Dir["../src_matlab_interface/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src_mex/" + File.basename(filename);
end
lst = Dir["../src_matlab_interface/*.hh"]
lst.each do |filename|
  FileUtils.cp filename, "./src_mex/" + File.basename(filename);
end

FileUtils.rm_rf   "tests"
FileUtils.mkdir_p "tests"
lst = Dir["../tests_matlab/*.m"]
lst.each do |filename|
  FileUtils.cp filename, "./tests/" + File.basename(filename);
end

FileUtils.cp "../license.txt", "license.txt"
