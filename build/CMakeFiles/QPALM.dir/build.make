# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ben/Documents/Projects/QPALM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ben/Documents/Projects/QPALM/build

# Include any dependencies generated for this target.
include CMakeFiles/QPALM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/QPALM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/QPALM.dir/flags.make

CMakeFiles/QPALM.dir/src/newton.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/newton.c.o: ../src/newton.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/QPALM.dir/src/newton.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/newton.c.o   -c /home/ben/Documents/Projects/QPALM/src/newton.c

CMakeFiles/QPALM.dir/src/newton.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/newton.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/newton.c > CMakeFiles/QPALM.dir/src/newton.c.i

CMakeFiles/QPALM.dir/src/newton.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/newton.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/newton.c -o CMakeFiles/QPALM.dir/src/newton.c.s

CMakeFiles/QPALM.dir/src/newton.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/newton.c.o.requires

CMakeFiles/QPALM.dir/src/newton.c.o.provides: CMakeFiles/QPALM.dir/src/newton.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/newton.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/newton.c.o.provides

CMakeFiles/QPALM.dir/src/newton.c.o.provides.build: CMakeFiles/QPALM.dir/src/newton.c.o


CMakeFiles/QPALM.dir/src/util.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/util.c.o: ../src/util.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/QPALM.dir/src/util.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/util.c.o   -c /home/ben/Documents/Projects/QPALM/src/util.c

CMakeFiles/QPALM.dir/src/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/util.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/util.c > CMakeFiles/QPALM.dir/src/util.c.i

CMakeFiles/QPALM.dir/src/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/util.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/util.c -o CMakeFiles/QPALM.dir/src/util.c.s

CMakeFiles/QPALM.dir/src/util.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/util.c.o.requires

CMakeFiles/QPALM.dir/src/util.c.o.provides: CMakeFiles/QPALM.dir/src/util.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/util.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/util.c.o.provides

CMakeFiles/QPALM.dir/src/util.c.o.provides.build: CMakeFiles/QPALM.dir/src/util.c.o


CMakeFiles/QPALM.dir/src/qpalm.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/qpalm.c.o: ../src/qpalm.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/QPALM.dir/src/qpalm.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/qpalm.c.o   -c /home/ben/Documents/Projects/QPALM/src/qpalm.c

CMakeFiles/QPALM.dir/src/qpalm.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/qpalm.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/qpalm.c > CMakeFiles/QPALM.dir/src/qpalm.c.i

CMakeFiles/QPALM.dir/src/qpalm.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/qpalm.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/qpalm.c -o CMakeFiles/QPALM.dir/src/qpalm.c.s

CMakeFiles/QPALM.dir/src/qpalm.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/qpalm.c.o.requires

CMakeFiles/QPALM.dir/src/qpalm.c.o.provides: CMakeFiles/QPALM.dir/src/qpalm.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/qpalm.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/qpalm.c.o.provides

CMakeFiles/QPALM.dir/src/qpalm.c.o.provides.build: CMakeFiles/QPALM.dir/src/qpalm.c.o


CMakeFiles/QPALM.dir/src/scaling.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/scaling.c.o: ../src/scaling.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/QPALM.dir/src/scaling.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/scaling.c.o   -c /home/ben/Documents/Projects/QPALM/src/scaling.c

CMakeFiles/QPALM.dir/src/scaling.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/scaling.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/scaling.c > CMakeFiles/QPALM.dir/src/scaling.c.i

CMakeFiles/QPALM.dir/src/scaling.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/scaling.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/scaling.c -o CMakeFiles/QPALM.dir/src/scaling.c.s

CMakeFiles/QPALM.dir/src/scaling.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/scaling.c.o.requires

CMakeFiles/QPALM.dir/src/scaling.c.o.provides: CMakeFiles/QPALM.dir/src/scaling.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/scaling.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/scaling.c.o.provides

CMakeFiles/QPALM.dir/src/scaling.c.o.provides.build: CMakeFiles/QPALM.dir/src/scaling.c.o


CMakeFiles/QPALM.dir/src/linesearch.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/linesearch.c.o: ../src/linesearch.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/QPALM.dir/src/linesearch.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/linesearch.c.o   -c /home/ben/Documents/Projects/QPALM/src/linesearch.c

CMakeFiles/QPALM.dir/src/linesearch.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/linesearch.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/linesearch.c > CMakeFiles/QPALM.dir/src/linesearch.c.i

CMakeFiles/QPALM.dir/src/linesearch.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/linesearch.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/linesearch.c -o CMakeFiles/QPALM.dir/src/linesearch.c.s

CMakeFiles/QPALM.dir/src/linesearch.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/linesearch.c.o.requires

CMakeFiles/QPALM.dir/src/linesearch.c.o.provides: CMakeFiles/QPALM.dir/src/linesearch.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/linesearch.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/linesearch.c.o.provides

CMakeFiles/QPALM.dir/src/linesearch.c.o.provides.build: CMakeFiles/QPALM.dir/src/linesearch.c.o


CMakeFiles/QPALM.dir/src/lin_alg.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/lin_alg.c.o: ../src/lin_alg.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/QPALM.dir/src/lin_alg.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/lin_alg.c.o   -c /home/ben/Documents/Projects/QPALM/src/lin_alg.c

CMakeFiles/QPALM.dir/src/lin_alg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/lin_alg.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/lin_alg.c > CMakeFiles/QPALM.dir/src/lin_alg.c.i

CMakeFiles/QPALM.dir/src/lin_alg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/lin_alg.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/lin_alg.c -o CMakeFiles/QPALM.dir/src/lin_alg.c.s

CMakeFiles/QPALM.dir/src/lin_alg.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/lin_alg.c.o.requires

CMakeFiles/QPALM.dir/src/lin_alg.c.o.provides: CMakeFiles/QPALM.dir/src/lin_alg.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/lin_alg.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/lin_alg.c.o.provides

CMakeFiles/QPALM.dir/src/lin_alg.c.o.provides.build: CMakeFiles/QPALM.dir/src/lin_alg.c.o


CMakeFiles/QPALM.dir/src/termination.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/termination.c.o: ../src/termination.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/QPALM.dir/src/termination.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/termination.c.o   -c /home/ben/Documents/Projects/QPALM/src/termination.c

CMakeFiles/QPALM.dir/src/termination.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/termination.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/termination.c > CMakeFiles/QPALM.dir/src/termination.c.i

CMakeFiles/QPALM.dir/src/termination.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/termination.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/termination.c -o CMakeFiles/QPALM.dir/src/termination.c.s

CMakeFiles/QPALM.dir/src/termination.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/termination.c.o.requires

CMakeFiles/QPALM.dir/src/termination.c.o.provides: CMakeFiles/QPALM.dir/src/termination.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/termination.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/termination.c.o.provides

CMakeFiles/QPALM.dir/src/termination.c.o.provides.build: CMakeFiles/QPALM.dir/src/termination.c.o


CMakeFiles/QPALM.dir/src/cholmod_interface.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/cholmod_interface.c.o: ../src/cholmod_interface.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/QPALM.dir/src/cholmod_interface.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/cholmod_interface.c.o   -c /home/ben/Documents/Projects/QPALM/src/cholmod_interface.c

CMakeFiles/QPALM.dir/src/cholmod_interface.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/cholmod_interface.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/cholmod_interface.c > CMakeFiles/QPALM.dir/src/cholmod_interface.c.i

CMakeFiles/QPALM.dir/src/cholmod_interface.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/cholmod_interface.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/cholmod_interface.c -o CMakeFiles/QPALM.dir/src/cholmod_interface.c.s

CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.requires

CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.provides: CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.provides

CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.provides.build: CMakeFiles/QPALM.dir/src/cholmod_interface.c.o


CMakeFiles/QPALM.dir/src/nonconvex.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/nonconvex.c.o: ../src/nonconvex.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/QPALM.dir/src/nonconvex.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/nonconvex.c.o   -c /home/ben/Documents/Projects/QPALM/src/nonconvex.c

CMakeFiles/QPALM.dir/src/nonconvex.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/nonconvex.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/nonconvex.c > CMakeFiles/QPALM.dir/src/nonconvex.c.i

CMakeFiles/QPALM.dir/src/nonconvex.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/nonconvex.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/nonconvex.c -o CMakeFiles/QPALM.dir/src/nonconvex.c.s

CMakeFiles/QPALM.dir/src/nonconvex.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/nonconvex.c.o.requires

CMakeFiles/QPALM.dir/src/nonconvex.c.o.provides: CMakeFiles/QPALM.dir/src/nonconvex.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/nonconvex.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/nonconvex.c.o.provides

CMakeFiles/QPALM.dir/src/nonconvex.c.o.provides.build: CMakeFiles/QPALM.dir/src/nonconvex.c.o


CMakeFiles/QPALM.dir/src/validate.c.o: CMakeFiles/QPALM.dir/flags.make
CMakeFiles/QPALM.dir/src/validate.c.o: ../src/validate.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/QPALM.dir/src/validate.c.o"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/QPALM.dir/src/validate.c.o   -c /home/ben/Documents/Projects/QPALM/src/validate.c

CMakeFiles/QPALM.dir/src/validate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QPALM.dir/src/validate.c.i"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ben/Documents/Projects/QPALM/src/validate.c > CMakeFiles/QPALM.dir/src/validate.c.i

CMakeFiles/QPALM.dir/src/validate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QPALM.dir/src/validate.c.s"
	/usr/bin/gcc-7 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ben/Documents/Projects/QPALM/src/validate.c -o CMakeFiles/QPALM.dir/src/validate.c.s

CMakeFiles/QPALM.dir/src/validate.c.o.requires:

.PHONY : CMakeFiles/QPALM.dir/src/validate.c.o.requires

CMakeFiles/QPALM.dir/src/validate.c.o.provides: CMakeFiles/QPALM.dir/src/validate.c.o.requires
	$(MAKE) -f CMakeFiles/QPALM.dir/build.make CMakeFiles/QPALM.dir/src/validate.c.o.provides.build
.PHONY : CMakeFiles/QPALM.dir/src/validate.c.o.provides

CMakeFiles/QPALM.dir/src/validate.c.o.provides.build: CMakeFiles/QPALM.dir/src/validate.c.o


# Object files for target QPALM
QPALM_OBJECTS = \
"CMakeFiles/QPALM.dir/src/newton.c.o" \
"CMakeFiles/QPALM.dir/src/util.c.o" \
"CMakeFiles/QPALM.dir/src/qpalm.c.o" \
"CMakeFiles/QPALM.dir/src/scaling.c.o" \
"CMakeFiles/QPALM.dir/src/linesearch.c.o" \
"CMakeFiles/QPALM.dir/src/lin_alg.c.o" \
"CMakeFiles/QPALM.dir/src/termination.c.o" \
"CMakeFiles/QPALM.dir/src/cholmod_interface.c.o" \
"CMakeFiles/QPALM.dir/src/nonconvex.c.o" \
"CMakeFiles/QPALM.dir/src/validate.c.o"

# External object files for target QPALM
QPALM_EXTERNAL_OBJECTS =

libQPALM.a: CMakeFiles/QPALM.dir/src/newton.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/util.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/qpalm.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/scaling.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/linesearch.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/lin_alg.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/termination.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/cholmod_interface.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/nonconvex.c.o
libQPALM.a: CMakeFiles/QPALM.dir/src/validate.c.o
libQPALM.a: CMakeFiles/QPALM.dir/build.make
libQPALM.a: CMakeFiles/QPALM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ben/Documents/Projects/QPALM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking C static library libQPALM.a"
	$(CMAKE_COMMAND) -P CMakeFiles/QPALM.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/QPALM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/QPALM.dir/build: libQPALM.a

.PHONY : CMakeFiles/QPALM.dir/build

CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/newton.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/util.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/qpalm.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/scaling.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/linesearch.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/lin_alg.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/termination.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/cholmod_interface.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/nonconvex.c.o.requires
CMakeFiles/QPALM.dir/requires: CMakeFiles/QPALM.dir/src/validate.c.o.requires

.PHONY : CMakeFiles/QPALM.dir/requires

CMakeFiles/QPALM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/QPALM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/QPALM.dir/clean

CMakeFiles/QPALM.dir/depend:
	cd /home/ben/Documents/Projects/QPALM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ben/Documents/Projects/QPALM /home/ben/Documents/Projects/QPALM /home/ben/Documents/Projects/QPALM/build /home/ben/Documents/Projects/QPALM/build /home/ben/Documents/Projects/QPALM/build/CMakeFiles/QPALM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/QPALM.dir/depend

