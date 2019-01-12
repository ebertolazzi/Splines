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

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLOBBER.include []
CLEAN.exclude('**/[cC][oO][rR][eE]')

task :default => [:build]

LIB_NAME="Splines"

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

desc "compile for UNIX/OSX [default GC='./GC']"
task :build, [:gc_dir] do |t, args|

  args.with_defaults( :gc_dir => './GC' )

  if args.gc_dir == './GC' then
    puts "\n\nBuild submodule GenericContainer".green
    FileUtils.rm_rf "GC"
    sh "git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC"
    FileUtils.cd "GC"
    FileUtils.cp "../CMakeLists-cflags.txt", "CMakeLists-cflags.txt"
    sh "rake build"
    FileUtils.cd ".."
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  FileUtils.rm_rf   "lib"
  FileUtils.mkdir_p "lib"
  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "\n\nPrepare #{LIB_NAME} project".green
  sh "cmake -DCMAKE_INSTALL_PREFIX:PATH=lib -DGC_DIR:PATH=\"#{args.gc_dir}\" .."

  puts "\n\nCompile #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug  --target install'
  FileList["*#{LIB_NAME}*.*"].each do |f|
    puts "Copying #{f}".yellow
    ext = File.extname(f);
    FileUtils.cp f, "../lib/#{File.basename(f,ext)}_debug#{ext}"
  end

  puts "\n\nCompile Splines Release".green
  sh 'cmake --build . --config Release --target install'
  FileList["*#{LIB_NAME}*.*"].each do |f|
    puts "Copying #{f}".yellow
    FileUtils.cp f, "../lib/#{File.basename(f)}"
  end

  puts "\n\nCopy include".green
  FileUtils.cp_r "lib/lib/include", "../lib/include"

  FileUtils.cd '..'

end

desc "compile for Visual Studio [default year=2017 bits=x64 gc_dir='./GC']"
task :build_win, [:year, :bits, :gc_dir] do |t, args|
  args.with_defaults( :year => "2017", :bits => "x64", :gc_dir => './GC' )

  if args.gc_dir == './GC' then
    puts "\n\nBuild submodule GenericContainer".green
    FileUtils.rm_rf "GC"
    sh "git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC"
    FileUtils.cd "GC"
    FileUtils.cp "../CMakeLists-cflags.txt", "CMakeLists-cflags.txt"
    sh "rake build_win[#{args.year},#{args.bits}]"
    FileUtils.cd ".."
  else
    puts "\n\nUse GenericContainer at #{args.gc_dir}".green
  end

  puts "\n\nPrepare #{LIB_NAME} project".green
  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   "lib"
  FileUtils.mkdir_p "lib"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " +
        ' -DCMAKE_INSTALL_PREFIX:PATH=lib' +
        " -DGC_DIR:PATH=\"#{args.gc_dir}\" .."

  win32_64 = ''
  case args.bits
  when /x64/
    win32_64 = ' Win64'
  end

  case args.year
  when "2010"
    sh 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
  when "2012"
    sh 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
  when "2013"
    sh 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
  when "2015"
    sh 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
  when "2017"
    sh 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
  else
    puts "Visual Studio year #{year} not supported!\n";
  end

  libname = "#{LIB_NAME}_vs#{args.year}_#{args.bits}"

  puts "\n\nCompile #{LIB_NAME} Debug".green
  sh 'cmake --build . --config Debug --target install'
  FileUtils.cp "Debug/#{LIB_NAME}.lib", "../lib/#{libname}_debug.lib"

  puts "\n\nCompile #{LIB_NAME} Release".green
  sh 'cmake --build . --config Release  --target install'
  FileUtils.cp "Release/#{LIB_NAME}.lib", "../lib/#{libname}.lib"  

  puts "\n\nCopy include".green
  FileUtils.cp_r "lib/lib/include", "../lib/include"

  FileUtils.cd '..'

end
