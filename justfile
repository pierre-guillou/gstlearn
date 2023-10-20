export CMAKE_EXPORT_COMPILE_COMMANDS := "ON"
export CMAKE_COLOR_DIAGNOSTICS := "ON"
export CMAKE_GENERATOR := "Ninja"
export CCACHE_DISABLE := "1"

bdir := "_build"

# apply clang-format on every source file
clang_format:
    fd -e cxx -e h -E sparselib -x clang-format -i

# check compiler warnings with clang-check
clang_check: (_pass_on_sources "clang-check")

# run clang-tidy on every source file
clang_tidy: (_pass_on_sources "clang-tidy")

# run Clang Static Analyzer (scan-build)
scan_build: (_build_dir "scan-build")
    scan-build -V cmake --build {{ bdir }}

# get a list of all sources from the compilation database with jq and
# apply a given clang tool on them
_pass_on_sources tool: _build_dir
    for f in $(jq -r '.[] | .file' < {{ bdir }}/compile_commands.json); do \
      just _clang_tool {{ tool }} $f; \
    done

# apply a given clang tool on a given file
_clang_tool tool file:
    {{ tool }} -p {{ bdir }} {{ file }}

# generate a brand new build directory
_build_dir wrapper="" cmake_opts="":
    [ -d {{ bdir }} ] && rm -rf {{ bdir }} || true
    mkdir {{ bdir }}
    {{ wrapper }} cmake -B {{ bdir }} {{ cmake_opts }}

# build with coverage support
build_coverage:
    CFLAGS="--coverage" CXXFLAGS="--coverage" just _build_dir "" "-DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=TRUE"
    cmake --build {{ bdir }} --parallel

# generate a coverage report
coverage_report:
    grcov . -s . --binary-path {{ bdir }} -t html --branch --ignore-not-existing -o {{ bdir }}/coverage
    xdg-open {{ bdir }}/coverage/index.html
