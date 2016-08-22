@echo off
@rmdir /S GC ; git clone --depth 1 git@github.com:ebertolazzi/GenericContainer.git GC
@pushd win_compile
@call sp-compile-vs2013
@call sp-compile-vs2015
@call sp-install-vs2013
@call sp-install-vs2015
@popd
