%w(colorize fileutils rake/clean).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end
require 'shellwords'

# avoid to remove file "core" (in Eigen inclusion)
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

#
# Check for a configuration file on a upper directory.
# This permits to use a unique configuration file for
# large projects.
# On a local project use the default in this file.
#
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the root of the local project
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__))) then
  # found in the upper project
  require_relative '../../Rakefile_configure.rb'
else
  #-------------------------
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
  #-------------------------
end

#    ___  ____
#   / _ \/ ___|
#  | | | \___ \
#  | |_| |___) |
#   \___/|____/
#
case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux|cygwin/ # cygwin compile as a linux system
  OS = :linux
when /msys/
  # msys2 envirorment to compile with MINGW
  OS = :mingw
else # assume windows
  OS = :win
end
def build_type
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

def install_prefix
  File.expand_path('lib', __dir__)
end

def project_root
  File.expand_path(__dir__)
end

def native_build_command(action)
  case OS
  when :mac, :linux, :mingw
    ['bash', File.expand_path('build.sh', __dir__), action, build_type, '-p', install_prefix]
  when :win
    ['pwsh', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', File.expand_path('build.ps1', __dir__), action, build_type, '-p', install_prefix]
  else
    raise "Unsupported platform #{OS}"
  end
end

desc "default task --> build"
task :default => :build

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end

#   ____  _   _ _   _
#  |  _ \| | | | \ | |
#  | |_) | | | |  \| |
#  |  _ <| |_| | |\  |
#  |_| \_\\___/|_| \_|
#
desc "build and run all tests"
task :run do
  puts "Run tests".green
  yellow_sh(*native_build_command('test'))
end

desc "run tests"
task :test do
  puts "Test".green
  yellow_sh(*native_build_command('test'))
end

#   ____  _   _ ___ _     ____
#  | __ )| | | |_ _| |   |  _ \
#  |  _ \| | | || || |   | | | |
#  | |_) | |_| || || |___| |_| |
#  |____/ \___/|___|_____|____/
#
desc "build"
task :build do
  puts "Build".green
  yellow_sh(*native_build_command('install'))
end

desc "clean"
task :clean do
  case OS
  when :mac
    puts "Clean (osx)".green
    Rake::Task[:clean_osx].invoke
  when :linux
    puts "Clean (linux)".green
    Rake::Task[:clean_linux].invoke
  when :win
    puts "Clean (windows)".green
    Rake::Task[:clean_win].invoke
  when :mingw
    puts "Clean (mingw)".green
    Rake::Task[:clean_mingw].invoke
  else
    raise "Unsupported platform #{OS}"
  end
end

desc "default task --> build"
task :default => :build

def yellow_sh(*cmd)
  puts cmd.map { |part| Shellwords.escape(part) }.join(' ').yellow
  Dir.chdir(project_root) { sh(*cmd) }
end

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }
CLOBBER.include []

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
  puts "Package".green
  yellow_sh(*native_build_command('package'))
end
