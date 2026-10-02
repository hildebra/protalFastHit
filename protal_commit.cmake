# The commit the source is at, for the binaries' --version (src/Utilities/BuildInfo.cpp): writes OUTPUT with
# PROTAL_GIT_COMMIT (empty without git or outside a checkout, e.g. a source tarball) and PROTAL_GIT_MODIFIED (1 when
# what the binaries are built from has uncommitted changes). Run at every build by the target protal_commit; OUTPUT
# is only rewritten when it changes, so nothing else is recompiled. build_gtdb_database.py compares the commit with
# its own checkout, so that a run does not go on for hours with binaries of an older source.
# cmake -DSOURCE_DIR=<checkout> -DOUTPUT=<header> -P protal_commit.cmake
set(commit "")
set(modified 0)
find_program(GIT_EXECUTABLE git)
if (GIT_EXECUTABLE)
	execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} rev-parse HEAD
			OUTPUT_VARIABLE head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE failed)
	if (NOT failed AND head MATCHES "^[0-9a-f]+$")
		set(commit ${head})
		# The paths build_gtdb_database.py compares (BUILD_SOURCES), untracked files in them included.
		execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SOURCE_DIR} status --porcelain --
				src lib CMakeLists.txt protal_config.h.in protal_commit.cmake
				OUTPUT_VARIABLE changes ERROR_QUIET)
		if (changes)
			set(modified 1)
		endif()
	endif()
endif()
file(WRITE ${OUTPUT}.new "#pragma once\n#define PROTAL_GIT_COMMIT \"${commit}\"\n#define PROTAL_GIT_MODIFIED ${modified}\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${OUTPUT}.new ${OUTPUT})
file(REMOVE ${OUTPUT}.new)
