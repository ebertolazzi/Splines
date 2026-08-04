# frozen_string_literal: true

%w[colorize fileutils].each do |lib|
  begin
    require lib
  rescue LoadError
    warn "Install the #{lib} gem:\n $ (sudo) gem install #{lib}"
    exit 1
  end
end

require 'rake/clean'
require 'etc'
require 'shellwords'

# Avoid removing files named "core" from vendored headers.
CLEAN.clear_exclude.exclude { |fn| fn.pathmap('%f').downcase == 'core' }

config_paths = [
  File.expand_path('../Rakefile_configure.rb', __dir__),
  File.expand_path('../../Rakefile_configure.rb', __dir__)
]
config_file = config_paths.find { |path| File.exist?(path) }

if config_file
  require config_file
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
end

OS = case RUBY_PLATFORM
     when /darwin/ then :mac
     when /linux|cygwin/ then :linux
     when /msys/ then :mingw
     else :win
     end

PROJECT_ROOT = File.expand_path(__dir__)
BUILD_DIR    = File.join(PROJECT_ROOT, 'build')
INSTALL_DIR  = File.join(PROJECT_ROOT, 'lib')
BIN_DIR      = File.join(PROJECT_ROOT, 'bin')
ALLOW_NETWORK_FETCH = ENV.fetch('NETWORK_FETCH', 'ON').match?(/\A(1|on|true|yes)\z/i)

CMAKE_BUILD_PARALLEL_ARGS = begin
  if OS == :win
    ['--parallel']
  else
    ['--parallel', Etc.nprocessors.to_s]
  end
end

def cmake_bool(value)
  value ? 'ON' : 'OFF'
end

def build_configuration
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

def command_string(*cmd)
  cmd.flatten.map { |part| Shellwords.escape(part.to_s) }.join(' ')
end

def yellow_sh(*cmd)
  puts command_string(*cmd).yellow
  sh(*cmd)
end

def configure_args(enable_tests: false)
  [
    'cmake',
    '-G', 'Ninja',
    '-S', PROJECT_ROOT,
    '-B', BUILD_DIR,
    "-DCMAKE_BUILD_TYPE:STRING=#{build_configuration}",
    "-DCMAKE_INSTALL_PREFIX:PATH=#{INSTALL_DIR}",
    '-DCMAKE_INSTALL_LIBDIR:PATH=lib',
    '-DCMAKE_INSTALL_INCLUDEDIR:PATH=include',
    "-DBUILD_SHARED_LIBS:BOOL=#{cmake_bool(COMPILE_DYNAMIC)}",
    "-DBUILD_TESTING:BOOL=#{cmake_bool(enable_tests)}",
    '-DSPLINES_BUILD_BENCHMARKS:BOOL=OFF',
    '-DSPLINES_INSTALL:BOOL=ON',
    '-DSPLINES_UPDATE_3RDPARTY:BOOL=OFF',
    '-DSPLINES_COLLECT_DEPENDENCIES:BOOL=OFF',
    '-DSPLINES_POPULATE_TOOLBOX:BOOL=OFF',
    "-DSPLINES_ALLOW_NETWORK_FETCH:BOOL=#{cmake_bool(ALLOW_NETWORK_FETCH)}",
    "-DGENERIC_CONTAINER_ALLOW_NETWORK_FETCH:BOOL=#{cmake_bool(ALLOW_NETWORK_FETCH)}",
    '-DUTILS_UPDATE_3RDPARTY:BOOL=OFF'
  ]
end

def cmake_build_args(target)
  [
    'cmake', '--build', BUILD_DIR,
    '--config', build_configuration,
    '--target', target,
    *CMAKE_BUILD_PARALLEL_ARGS
  ]
end

def configure(enable_tests: false)
  FileUtils.mkdir_p(BUILD_DIR)
  yellow_sh(*configure_args(enable_tests: enable_tests))
end

desc 'default task: build library only'
task default: :build

desc 'reset git repository and remove ignored files'
task :git_clean do
  sh 'git reset --hard'
  sh 'git clean -d -x -f'
end

desc 'build and install the Splines library only, without tests or benchmarks'
task :build do
  puts "Build library only (#{OS})".green
  configure(enable_tests: false)
  yellow_sh(*cmake_build_args('install'))
end

desc 'build and run Splines tests'
task :run do
  puts 'Build and run tests'.green
  configure(enable_tests: true)
  yellow_sh(*cmake_build_args('Splines_all_tests'))
  yellow_sh('ctest', '--test-dir', BUILD_DIR, '--build-config', build_configuration, '--output-on-failure')
end

desc 'alias for run'
task test: :run

desc 'remove generated files'
task :clean do
  FileUtils.rm_rf(BUILD_DIR)
  FileUtils.rm_rf(INSTALL_DIR)
  FileUtils.rm_rf(BIN_DIR)
  FileUtils.rm_rf(File.join(PROJECT_ROOT, 'lib3rd'))
end

%i[osx linux mingw win].each do |platform|
  task "clean_#{platform}" => :clean
  task "build_#{platform}" => :build
end

desc 'build source/binary packages with CPack'
task :cpack do
  configure(enable_tests: false) unless Dir.exist?(BUILD_DIR)
  yellow_sh('cmake', '--build', BUILD_DIR, '--config', build_configuration, '--target', 'package')
end

CLEAN.include './**/*.o', './**/*.obj', './bin/**/example*', './build'
CLOBBER.include
