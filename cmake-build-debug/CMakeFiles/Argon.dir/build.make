# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Administrator\Desktop\KMS1\Argon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Argon.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Argon.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Argon.dir/flags.make

CMakeFiles/Argon.dir/argon.cpp.obj: CMakeFiles/Argon.dir/flags.make
CMakeFiles/Argon.dir/argon.cpp.obj: ../argon.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Argon.dir/argon.cpp.obj"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Argon.dir\argon.cpp.obj -c C:\Users\Administrator\Desktop\KMS1\Argon\argon.cpp

CMakeFiles/Argon.dir/argon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Argon.dir/argon.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Administrator\Desktop\KMS1\Argon\argon.cpp > CMakeFiles\Argon.dir\argon.cpp.i

CMakeFiles/Argon.dir/argon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Argon.dir/argon.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Administrator\Desktop\KMS1\Argon\argon.cpp -o CMakeFiles\Argon.dir\argon.cpp.s

# Object files for target Argon
Argon_OBJECTS = \
"CMakeFiles/Argon.dir/argon.cpp.obj"

# External object files for target Argon
Argon_EXTERNAL_OBJECTS =

Argon.exe: CMakeFiles/Argon.dir/argon.cpp.obj
Argon.exe: CMakeFiles/Argon.dir/build.make
Argon.exe: CMakeFiles/Argon.dir/linklibs.rsp
Argon.exe: CMakeFiles/Argon.dir/objects1.rsp
Argon.exe: CMakeFiles/Argon.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Argon.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Argon.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Argon.dir/build: Argon.exe

.PHONY : CMakeFiles/Argon.dir/build

CMakeFiles/Argon.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Argon.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Argon.dir/clean

CMakeFiles/Argon.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Administrator\Desktop\KMS1\Argon C:\Users\Administrator\Desktop\KMS1\Argon C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug C:\Users\Administrator\Desktop\KMS1\Argon\cmake-build-debug\CMakeFiles\Argon.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Argon.dir/depend
