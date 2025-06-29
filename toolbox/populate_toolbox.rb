#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf "src"
Dir.glob("bin/*.mex*").each { |file| File.delete(file)}

FileUtils.cp_r  "../src/.", "./src";
FileUtils.cp_r  "../cmake_utils/.",   "./cmake_utils";
FileUtils.rm_rf "./cmake_utils/.git"; # remove git struture
FileUtils.cp_r  "../submodules/quarticRootsFlocke/src/.",   "./src";
FileUtils.cp_r  "../submodules/UtilsLite/src/.",            "./src";
FileUtils.cp_r  "../submodules/GenericContainer/src/.",     "./src";
FileUtils.cp_r  "../submodules/GenericContainer/include/.", "./src";
FileUtils.cp    "../submodules/GenericContainer/src_matlab_interface/GenericContainerInterface_matlab.cc", "./src";

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

FileUtils.cp "../submodules/GenericContainer/src_matlab_interface/GenericContainerInterface_matlab.cc",
             "./src_mex/GenericContainerInterface_matlab.cc"

FileUtils.cp "../submodules/GenericContainer/src_json_interface/to_json.cc",
             "./src/to_json.cc"
FileUtils.cp "../submodules/GenericContainer/src_json_interface/from_json.cc",
             "./src/from_json.cc"
FileUtils.cp "../submodules/GenericContainer/src_toml_interface/to_toml.cc",
             "./src/to_toml.cc"
FileUtils.cp "../submodules/GenericContainer/src_toml_interface/from_toml.cc",
             "./src/from_toml.cc"

FileUtils.cp "../submodules/GenericContainer/src_yaml_interface/to_yaml.cc",
             "./src/to_yaml.cc"
FileUtils.cp "../submodules/GenericContainer/src_yaml_interface/from_yaml.cc",
             "./src/from_yaml.cc"

FileUtils.cp_r  "../submodules/GenericContainer/src_toml_interface/toml++/.", "./src/toml++";
FileUtils.cp_r  "../submodules/GenericContainer/src_json_interface/rapidjson/.", "./src/rapidjson";
FileUtils.cp_r  "../submodules/GenericContainer/src_yaml_interface/fkYAML/.", "./src/fkYAML";

FileUtils.cp "../license.txt", "license.txt"
