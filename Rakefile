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

desc "compile for Visual Studio [default year=2017, bits=x64, GC='./GC']"
task :build_GC, [:cmd] do |t, args|
  args.with_defaults( :cmd => "build_osx" )
  puts "\n\nBuild submodule GenericContainer".green
  FileUtils.rm_rf "GC"
  sh "git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC"
  FileUtils.cd "GC"
  FileUtils.cp "../CMakeLists-cflags.txt", "CMakeLists-cflags.txt"
  sh "rake #{args.cmd}"
  FileUtils.cd ".."
end

desc "compile for Visual Studio [default year=2017, bits=x64, GC='./GC']"
task :build_win, [:year, :bits, :gc_dir ] do |t, args|
  args.with_defaults( :year => "2017", :bits => "x64", :gc_dir => "./GC" )

  if args.gc_dir == './GC' then
    Rake::Task[:build_GC].invoke("build_win[#{args.year},#{args.bits}]")
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmake_cmd = win_vs(args.bits,args.year)
  cmake_cmd += " -DGC_DIR:VAR=#{args.gc_dir} "
  if COMPILE_EXECUTABLE then
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=false '
  end
  if COMPILE_DYNAMIC then
    cmake_cmd += ' -DBUILD_SHARED:VAR=true '
  else
    cmake_cmd += ' -DBUILD_SHARED:VAR=false '
  end

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  if COMPILE_DEBUG then
    sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end

  sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake  --build . --config Release  --target install '+PARALLEL
  FileUtils.cd '..'

end

desc "compile for OSX [default GC='./GC']"
task :build, [:gc_dir,:os] do |t, args|
  args.with_defaults( :gc_dir => "./GC" )

  if args.gc_dir == './GC' then
    Rake::Task[:build_GC].invoke("build_#{args.os}")
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmake_cmd = 'cmake '
  cmake_cmd += ' -DGC_DIR:VAR=#{args.gc_dir} '

  if COMPILE_EXECUTABLE then
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=false '
  end
  if COMPILE_DYNAMIC then
    cmake_cmd += '-DBUILD_SHARED:VAR=true '
  else
    cmake_cmd += '-DBUILD_SHARED:VAR=false '
  end

  if COMPILE_DEBUG then
    sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end
  sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --config Release --target install '+PARALLEL

  FileUtils.cd '..'
end

desc "compile for LINUX [default GC='./GC']"
task :build_linux, [:gc_dir ] do |t, args|
  args.with_defaults( :gc_dir => "./GC" )
  Rake::Task[:build].invoke(args.gc_dir,"linux")
end

desc "compile for OSX [default GC='./GC']"
task :build_osx, [:gc_dir ] do |t, args|
  args.with_defaults( :gc_dir => "./GC" )
  Rake::Task[:build].invoke(args.gc_dir,"osx")
end
