if File.exist?(File.expand_path('./cmake_utils/Rakefile_common.rb', File.dirname(__FILE__))) then
  require_relative "./cmake_utils/Rakefile_common.rb"
else
  require_relative "../Rakefile_common.rb"
end

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }
CLOBBER.include []

desc "compile for Visual Studio"
task :build_win do
  # check architecture
  case `where cl.exe`.chop
  when /(x64|amd64)\\cl\.exe/
    VS_ARCH = 'x64'
  when /(bin|x86|amd32)\\cl\.exe/
    VS_ARCH = 'x86'
  else
    raise RuntimeError, "Cannot determine architecture for Visual Studio".red
  end

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for SPLINES".yellow
  sh "cmake -G Ninja -DBITS:VAR=#{VS_ARCH} " + cmd_cmake_build() + ' ..'

  puts "compile with CMAKE for SPLINES".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL
  end

  FileUtils.cd '..'
end


desc "compile for OSX/LINUX/MINGW"
task :build_common do

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for SPLINES".yellow
  sh "cmake -G Ninja " + cmd_cmake_build() + ' ..'

  puts "compile with CMAKE for SPLINES".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release  --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

task :build_linux => :build_common do end
task :build_osx   => :build_common do end
task :build_mingw => :build_common do end

task :clean_common do
  FileUtils.rm_rf 'build'
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_osx   => :clean_common do end
task :clean_linux => :clean_common do end
task :clean_mingw => :clean_common do end
task :clean_win   => :clean_common do end

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for SPLINES".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
