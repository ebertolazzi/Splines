#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf "src"
Dir.glob("bin/*.mex*").each { |file| File.delete(file)}

FileUtils.cp_r  "../src/.", "./src";
FileUtils.cp_r  "../cmake_utils/.",   "./cmake_utils";
FileUtils.rm_rf "./cmake_utils/.git"; # remove git struture
FileUtils.cp_r  "../submodules/quarticRootsFlocke/src/.",   "./src";
FileUtils.cp_r  "../submodules/Utils/src/.",                "./src";
FileUtils.cp_r  "../submodules/GenericContainer/src/.",     "./src";
FileUtils.cp_r  "../submodules/GenericContainer/include/.", "./src";
# elimino dipendenze da Eigen
FileUtils.rm_rf "./src/Eigen";
FileUtils.rm_rf "./src/Utils_Poly.cc";
FileUtils.rm_rf "./src/Utils_GG2D.cc";
FileUtils.rm_rf "./src/Utils_HJPatternSearch.cc";
FileUtils.rm_rf "./src/Utils_NelderMead.cc";

#lst = Dir["../doc/*"]
#lst.each do |filename|
#  FileUtils.cp filename, "./doc/" + File.basename(filename);
#end

FileUtils.cp "../license.txt", "license.txt"

FileUtils.cp "../submodules/GenericContainer/src_matlab_interface/GenericContainerMatlabInterface.cc",
             "./src_mex/GenericContainerMatlabInterface.cc"

FileUtils.cp "../license.txt", "license.txt"
