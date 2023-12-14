

###########################################################
#       This file is just for convenience.
#       Use cmake .. in the build directory for compilation
###########################################################
#Tarball of Everything
tarball:
	tar cvzf SPINAS.tgz source/*.cpp include/* Makefile README LICENSE CMakeLists.txt tests/*cpp SM/*cpp SM/include/* SM/CH/main_22.c SM/CH/make_main SM/CH/models/*1.mdl scans/*.py scans/*.cpp scans/*.h


filestructure:
	tree . -o file_structure.txt -I 'build|scans|other|CH|.*'
