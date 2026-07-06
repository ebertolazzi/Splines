#!/usr/bin/env pwsh

# This script is used to build the project using CMake.
# Usage:
#   ./build.ps1 [option] [build_type] [-p install_prefix]
# Where:
#   option: "configure", "build", "test", "install", "package" (default: build)
#   build_type: "Debug", "Release", "RelWithDebInfo", "MinSizeRel" (default: Release)
#   install_prefix: path to install the project (default: /usr/local)
#
# Build directory is "build" by default. You can change it by setting the BUILD_DIR environment variable.
# Generator can be specified by setting the CMAKE_GENERATOR environment variable (default: Ninja).

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Show-Usage {
  @'
Usage:
  ./build.ps1 [option] [build_type] [-p install_prefix]

Where:
  option:         configure, build, test, install, package (default: build)
  build_type:     Debug, Release, RelWithDebInfo, MinSizeRel (default: Release)
  install_prefix: path to install the project (default: /usr/local)

Environment:
  BUILD_DIR        Build directory (default: build)
  CMAKE_GENERATOR  CMake generator (default: Ninja)

Examples:
  ./build.ps1 configure Debug
  ./build.ps1 build Release
  ./build.ps1 test
  ./build.ps1 install Release -p "$HOME/.local"
  ./build.ps1 package
'@
}

function Stop-WithUsage {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Message
  )

  [Console]::Error.WriteLine("build.ps1: $Message")
  [Console]::Error.WriteLine("")
  [Console]::Error.WriteLine((Show-Usage))
  exit 1
}

function Invoke-CheckedCommand {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Command,

    [Parameter(Mandatory = $true)]
    [string[]]$CommandArgs
  )

  & $Command @CommandArgs
  if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
  }
}

$ProjectRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
  $ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
}
Set-Location $ProjectRoot

$Option = "build"
$BuildType = "Release"
$InstallPrefix = "/usr/local"
$OptionSeen = $false
$BuildTypeSeen = $false

$ValidOptions = @("configure", "build", "test", "install", "package")
$ValidBuildTypes = @("Debug", "Release", "RelWithDebInfo", "MinSizeRel")

for ($I = 0; $I -lt $args.Count; $I++) {
  $Arg = $args[$I]

  switch ($Arg) {
    { $_ -cin @("-h", "--help") } {
      Show-Usage
      exit 0
    }
    { $_ -cin @("-p", "--prefix") } {
      if ($I + 1 -ge $args.Count) {
        Stop-WithUsage "missing argument for $Arg"
      }
      $InstallPrefix = $args[$I + 1]
      $I++
      continue
    }
    { $_ -cin $ValidOptions } {
      if ($OptionSeen) {
        Stop-WithUsage "option specified more than once"
      }
      $Option = $Arg
      $OptionSeen = $true
      continue
    }
    { $_ -cin $ValidBuildTypes } {
      if ($BuildTypeSeen) {
        Stop-WithUsage "build_type specified more than once"
      }
      $BuildType = $Arg
      $BuildTypeSeen = $true
      continue
    }
    default {
      Stop-WithUsage "unknown argument: $Arg"
    }
  }
}

$BuildDir = $env:BUILD_DIR
if ([string]::IsNullOrWhiteSpace($BuildDir)) {
  $BuildDir = "build"
}

$Generator = $env:CMAKE_GENERATOR
if ([string]::IsNullOrWhiteSpace($Generator)) {
  $Generator = "Ninja"
}

function Invoke-Configure {
  Invoke-CheckedCommand "cmake" @(
    "-S", ".",
    "-B", $BuildDir,
    "-G", $Generator,
    "-DCMAKE_BUILD_TYPE=$BuildType",
    "-DCMAKE_INSTALL_PREFIX=$InstallPrefix"
  )
}

function Ensure-Configured {
  $CachePath = Join-Path $BuildDir "CMakeCache.txt"
  if (-not (Test-Path -LiteralPath $CachePath -PathType Leaf)) {
    Invoke-Configure
  }
}

switch ($Option) {
  "configure" {
    Invoke-Configure
  }
  "build" {
    Ensure-Configured
    Invoke-CheckedCommand "cmake" @("--build", $BuildDir, "--config", $BuildType)
  }
  "test" {
    Ensure-Configured
    Invoke-CheckedCommand "cmake" @("--build", $BuildDir, "--config", $BuildType)
    Invoke-CheckedCommand "ctest" @("--test-dir", $BuildDir, "--build-config", $BuildType)
  }
  "install" {
    Ensure-Configured
    Invoke-CheckedCommand "cmake" @("--build", $BuildDir, "--config", $BuildType)
    Invoke-CheckedCommand "cmake" @("--install", $BuildDir, "--config", $BuildType, "--prefix", $InstallPrefix)
  }
  "package" {
    Ensure-Configured
    Invoke-CheckedCommand "cmake" @("--build", $BuildDir, "--config", $BuildType, "--target", "package")
  }
}
