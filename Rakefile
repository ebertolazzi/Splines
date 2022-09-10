%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require 'rake/clean'

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build", "./lib", "./lib3rd"]
CLOBBER.include []
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux/
  OS = :linux
when /cygwin|mswin|mingw|bccwin|wince|emx/
  OS = :win
when /msys/
  OS = :win
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=ON '
else
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=OFF '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=ON '
else
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=OFF '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end

desc "default task --> build"
task :default => :build

FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/Utils/cmake/CMakeLists-cflags.txt'
FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/quarticRootsFlocke/cmake/CMakeLists-cflags.txt'
FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/GenericContainer/cmake/CMakeLists-cflags.txt'

desc "run tests"
task :test do
  FileUtils.cd "build"
  sh 'ctest --output-on-failure'
  FileUtils.cd '..'
end

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

desc "build SPLINES"
task :build do
  case OS
  when :mac
    puts "SPLINES build (osx)".green
    Rake::Task[:build_osx].invoke
  when :linux
    puts "SPLINES build (linux)".green
    Rake::Task[:build_linux].invoke
  when :win
    puts "SPLINES build (windows)".green
    Rake::Task[:build_win].invoke
  end
end

desc "clean SPLINES"
task :clean do
  case OS
  when :mac
    puts "SPLINES clean (osx)".green
    Rake::Task[:clean_osx].invoke
  when :linux
    puts "SPLINES clean (linux)".green
    Rake::Task[:clean_linux].invoke
  when :win
    puts "SPLINES clean (windows)".green
    Rake::Task[:clean_win].invoke
  end
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMAKE for SPLINES".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for SPLINES".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for SPLINES".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end


desc "compile for OSX/LINUX"
task :build_common do

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
    sh 'cmake --build . --config Release  --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for SPLINES".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end

desc "compile for LINUX"
task :build_linux do
  Rake::Task[:build_common].invoke()
end

desc "compile for OSX"
task :build_osx do
  Rake::Task[:build_common].invoke()
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
