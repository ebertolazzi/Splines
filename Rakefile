#
#
#

%w( colorize rake fileutils ).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require "rake/clean"
require_relative "./Rakefile_common.rb"

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLOBBER.include []
CLEAN.exclude('**/[cC][oO][rR][eE]')

file_base = File.expand_path(File.dirname(__FILE__)).to_s+'/lib'

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DBUILD_EXECUTABLE:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_EXECUTABLE:VAR=false '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DBUILD_SHARED:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_SHARED:VAR=false '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=WARNING '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=WARNING '
end
cmd_cmake_build += " -DINSTALL_HERE:VAR=true "

FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/Utils/CMakeLists-cflags.txt'
FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/CMakeLists-cflags.txt'
FileUtils.cp 'CMakeLists-cflags.txt', 'submodules/GenericContainer/CMakeLists-cflags.txt'

desc "run tests"
task :run do
  sh "./bin/test1"
  sh "./bin/test2"
  sh "./bin/test3"
  sh "./bin/test4"
  sh "./bin/test5"
  sh "./bin/test6"
  ##sh "./bin/test7"
  sh "./bin/test8"
  sh "./bin/test9"
end

desc "run tests"
task :run_win do
  sh "./bin/Release/test1"
  sh "./bin/Release/test2"
  sh "./bin/Release/test3"
  sh "./bin/Release/test4"
  sh "./bin/Release/test5"
  sh "./bin/Release/test6"
  ##sh "./bin/Release/test7"
  sh "./bin/Release/test8"
  sh "./bin/Release/test9"
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  #cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build
  cmd_cmake = 'cmake -G "NMake Makefiles" ' + cmd_cmake_build

  puts "run CMAKE for SPLINES".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for SPLINES".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end


desc "compile for OSX"
task :build, [:os] do |t, args|

  args.with_defaults( :os => "osx" )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build

  puts "run CMAKE for SPLINES".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for SPLINES".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for LINUX"
task :build_linux do
  Rake::Task[:build].invoke("linux")
end

desc "compile for OSX"
task :build_osx do
  Rake::Task[:build].invoke("osx")
end

task :clean_osx do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_linux do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end
