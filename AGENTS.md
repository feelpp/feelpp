# Feel++ Coding & Tooling Primer

This primer keeps humans and LLMs aligned with the Feel++ coding rules. Follow it for any change unless a
more specific directory guide overrides it.

## Build & Test
- Configure with `cmake --preset default`; build via `cmake --build build/default -j`.
- Export `CMAKE_EXPORT_COMPILE_COMMANDS=ON` in presets or `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON` to enable clang-tidy.
- Run focused tests first (`ctest -R <target>` or `pytest` in `python/pyfeelpp`). Avoid rebuilding everything unless needed.

## Formatting
- Use `.clang-format` (Allman braces, 4-space indent, no hard column limit, pointers on the type).
- Never reformat third_party/ or external/ trees.
- Run `clang-format` only on modified files: `clang-format -i path/to/file.cpp`.
- Namespaces stay flush-left and use the compact C++17 form: `namespace feelpp::mesh {`.

## Naming (clang-tidy enforced)
- Namespaces: `lower_snake`
- Classes/structs/enums/concepts: `PascalCase`
- Functions/methods: `PascalCase`
- Members (static or non-static): prefer `M_PascalCase`; legacy `PascalCase_` is also accepted.
- Variables (locals, parameters): `lower_case`
- Constants/macros: `UPPER_CASE`

> Namespaces flush-left + compact, member `M_PascalCase` or `PascalCase_`, no mass reformat.

## C++20 Essentials
- Prefer `constexpr`, `consteval`, `[[nodiscard]]`, `noexcept` when they express intent.
- Use `std::span`/`std::string_view` instead of raw pointers where ownership stays external.
- Concepts > SFINAE for new APIs; mark deleted/defaulted special members explicitly.
- Use `std::filesystem` (not Boost) for paths.
- Do **not** `using namespace std;`.

## Exceptions & RAII
- Throw `std::runtime_error` or project-specific derived classes with clear messages.
- Use RAII for resources (MPI communicators, Kokkos views, file handles). Avoid naked `new/delete`.
- Prefer `std::unique_ptr`/`std::shared_ptr` when ownership is dynamic.

## Logging & Diagnostics
- Use Feel++ logging helpers/macros; avoid raw `std::cout` in libraries.
- Keep debug output behind verbosity guards. Never print secrets or tokens.

## HPC Performance Guardrails
- No heap allocation, file I/O, or blocking synchronization inside hot kernels (`KOKKOS_LAMBDA`, OpenMP loops, MPI collectives).
- Minimize MPI global barriers; state them explicitly if unavoidable.
- Capture the smallest set of values in Kokkos lambdas.
- Prefer `const`/`span` views to signal read-only data.

## Git Hygiene
- Keep diffs focused; do not refactor unrelated files.
- Run pre-commit (`pre-commit run --all-files`) before pushing.
- Commit messages: `component: concise summary`.
- PRs must pass clang-format, clang-tidy, codespell, and relevant tests.

## Good / Bad Examples
````cpp
// ✅ Good
namespace feelpp::mesh
{
class MeshRefiner
{
  public:
    void Refine(Mesh const & mesh, int levels);

  private:
    double M_TargetRatio{0.25};
};
}

// ❌ Bad
using namespace std;
namespace FeelPP {
class mesh_refiner {
  double _targetratio;
  void refine(mesh const& mesh,int levels){ cout << "refine" << endl; }
};
}
````
````cpp
// ✅ Kokkos
Kokkos::parallel_for("update", range, KOKKOS_LAMBDA(int i) {
    state(i) = alpha * input(i);
});

// ❌ Kokkos
Kokkos::parallel_for("update", range, KOKKOS_LAMBDA(int i) {
    std::ofstream file("out.txt");
    file << input(i);
});
````

