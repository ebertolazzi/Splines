#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src/Utils/fmt"
FileUtils.mkdir_p "src/Utils/zstream"

lst = Dir["../src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../src/Splines/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Splines/" + File.basename(filename);
end


lst = Dir["../submodules/GenericContainer/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/GenericContainer/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/Utils/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/Utils/fmt/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/fmt/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/Utils/zstream/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/zstream/" + File.basename(filename);
end

lst = Dir["../submodules/quarticRootsFlocke/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

lst = Dir["../submodules/quarticRootsFlocke/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end

FileUtils.cp "../submodules/GenericContainer/src_matlab_interface/GenericContainerMatlabInterface.cc",
             "./src/GenericContainerMatlabInterface.cc"
FileUtils.cp "../submodules/GenericContainer/src_matlab_interface/GenericContainerMatlabInterface.hh",
             "./src/GenericContainerMatlabInterface.hh"

FileUtils.cp "../license.txt", "license.txt"
