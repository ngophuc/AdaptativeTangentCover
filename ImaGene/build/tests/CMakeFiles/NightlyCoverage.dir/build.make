# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build

# Utility rule file for NightlyCoverage.

# Include the progress variables for this target.
include tests/CMakeFiles/NightlyCoverage.dir/progress.make

tests/CMakeFiles/NightlyCoverage:
	cd /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build/tests && /usr/local/Cellar/cmake/3.19.2/bin/ctest -D NightlyCoverage

NightlyCoverage: tests/CMakeFiles/NightlyCoverage
NightlyCoverage: tests/CMakeFiles/NightlyCoverage.dir/build.make

.PHONY : NightlyCoverage

# Rule to build all files generated by this target.
tests/CMakeFiles/NightlyCoverage.dir/build: NightlyCoverage

.PHONY : tests/CMakeFiles/NightlyCoverage.dir/build

tests/CMakeFiles/NightlyCoverage.dir/clean:
	cd /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/NightlyCoverage.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/NightlyCoverage.dir/clean

tests/CMakeFiles/NightlyCoverage.dir/depend:
	cd /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/tests /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build/tests /Users/hngo/Downloads/AdaptiveTangentCover/ImaGene/build/tests/CMakeFiles/NightlyCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/NightlyCoverage.dir/depend

