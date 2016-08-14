@echo off
@pushd win_compile
@call sp-compile-vs2013
@call sp-compile-vs2015
@call sp-install-vs2013
@call sp-install-vs2015
@popd
