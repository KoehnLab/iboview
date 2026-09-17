# AGENTS.md

Notes for coding agents (and human contributors) working on IboView.

## What this project is

IboView is a molecular-orbital viewer built around Intrinsic Bonding Orbitals
(IAO/IBO) analysis. It lets chemists load wave functions from external
quantum-chemistry programs (Molpro, Orca, Turbomole, ...) or compute simple
ones internally, then interactively visualize orbitals, bond orders, and
oxidation-state analysis in a 3D OpenGL view. It has its own JavaScript-based
scripting/automation layer for driving the UI and saving/restoring view state.

Upstream is unmaintained; this repo is a fork that keeps it building against
current toolchains (currently Qt6/C++20) and fixes issues as they come up.

## Build

Qt6 + qmake, C++20. See `README.md` for per-distro package lists.

```bash
mkdir build && cd build
qmake6 ../main.pro   # or: qmake, if it already points at Qt6
make -j$(nproc)
```

Boost and an external BLAS/LAPACK (e.g. MKL) are optional; without them the
build falls back to a bundled Eigen-based emulation (slower, but functional
for viewer-only use — see comments at the top of `main.pro`).

There is no display/GPU in most agent sandboxes. A clean compile+link is
verifiable there; anything about actual rendering, window behavior, or
interactive UI needs a human to confirm on a real machine.

## Source layout

- `src/IboView/` — the application: main window (`IvMain.*`), document/data
  model (`IvDocument.*`, `IvDataSet.*`), the 3D view and low-level GL wrappers
  (`IvView3D.*`, `IvGl.*`), the scripting engine (`IvScript.*`), volume/orbital
  rendering caches (`IvVolumeDataSet.*`), and the various dialog forms.
- `src/MicroScf/` — the embedded, self-contained SCF/DFT engine (basis sets,
  integrals driver, Fock builds, IAO construction, orbital localization).
- `src/IrCore/` — integral evaluation routines (ERIs, ECPs) used both by
  `MicroScf` and directly by `IvIao.cpp`/`IvAnalysis`-related code for
  post-processing imported wave functions (IAO/IBO/bond-order/EOS analysis
  need integrals over the actual basis even when the wave function itself
  was computed elsewhere).
- `src/Common/` — shared low-level utilities (linear algebra, I/O, grids,
  DIIS, atom data, ...), prefixed `Cx*`.
- `src/eigen/`, `src/pugixml/`, `src/GL/` (bundled GLEW), `src/External/` —
  vendored third-party code. Don't "clean these up" or modernize them as a
  drive-by; they're intentionally kept close to upstream.
- `resources/`, `example-data/`, `example-scripts/` — Qt resources, sample
  wave function files, and sample `.js` state/preset scripts.

### Generated property boilerplate

Files named `prop_*.h.inl` / `prop_*.cpp.inl` (e.g. `prop_FView3d.h.inl`) are
generated (a `make_properties.py`-style tool, not present in-tree) from a
simpler spec and `#include`d into their owning class. They carry a literal
`// this code is GENERATED. Do not change it directly!` banner — respect it;
if a property's behavior needs to change, find the underlying generator/spec
or make the smallest possible hand-patch and say so, don't silently diverge
from what the generator would produce.

## Known pitfalls (found the hard way — read before touching the render path)

- **`QOpenGLWidget`'s GL context is only guaranteed current inside
  `initializeGL()`/`paintGL()`/`resizeGL()`.** Any code reachable from a
  plain Qt slot or event handler (menu actions, mouse handlers, property
  setters, etc.) that touches GL objects directly needs an explicit
  `makeCurrent()` first. This app used to run on Qt5's `QGLWidget`, which
  tolerated GL calls from almost anywhere; that tolerance is gone.
- **But `QOpenGLWidget::makeCurrent()` also silently rebinds the widget's own
  internal backing framebuffer**, not just the GL context. Calling it from a
  function that's *also* reached mid-frame from deep inside the normal render
  pipeline (where a specific custom FBO — `pMainFbo`, `pDpLayer0`,
  `pPickFbo` — is already deliberately bound) will yank that binding away and
  send whatever draws next into the wrong, invisible target for that frame.
  If a helper function is called both from outside `paintGL` (needs
  `makeCurrent()`) and from inside it (must not touch context/FBO state),
  put the `makeCurrent()` call at the outside-`paintGL` call site, not inside
  the shared helper. `FViewImpl::RenderPickBuffer()` shows the safe pattern:
  `makeCurrent()` immediately followed by rebinding the FBO it actually needs.
- **`FView3d::paintGL()` has a `QMutex::tryLock()`-based re-entrancy guard**
  (its own long-standing comment calls it "a horribly ugly hack"). A burst of
  property changes (e.g. a script or preset setting many `Q_PROPERTY`s in a
  row, each wired to `update()`) can cause `paintGL()` to re-enter itself; a
  failed `tryLock()` now reschedules via `update()` rather than silently
  dropping the frame — if you touch this function, keep that property.
- **The scripting engine is `QJSEngine`** (Qt6; `QtScript`/`QScriptEngine` is
  gone). It has no way to bind a raw C++ callback as a variadic script
  function — see `FScriptGlobals` in `IvScript.cpp` for the
  `Q_INVOKABLE`-method + JS-wrapper-function pattern used for `format`,
  `join`, `Vec3`, `FreeLine`. It also has no live-prototype mechanism for a
  custom C++ value type — `FVec3d` crossing into script is a plain JS object
  built fresh on every crossing, not a shared, connected instance.
- **Qt6 module list**: `openglwidgets`, `qml` (for `QJSEngine`), and
  `svgwidgets` are separate modules from `opengl`/`script`/`svg` in Qt5 —
  see `main.pro`.
- **`QFlags` no longer converts implicitly from a plain `0`.** A default
  argument like `Qt::WindowFlags flags = 0` won't compile under Qt6; use the
  flag type's own default constructor (`Qt::WindowFlags()`).

## Debugging segfaults and memory corruption

- **A crash's backtrace shows where corruption was *detected*, not
  necessarily where it happened.** Heap corruption (double-free,
  use-after-free, buffer overflow) can crash much later than the code that
  caused it — an unrelated-looking destructor, allocator call, or even a
  perfectly ordinary `delete` can be where a previously-corrupted object
  finally gets touched. Treat the function named in a plain debugger
  backtrace as the symptom site, not necessarily the cause, especially if a
  fix there doesn't change the crash at all.
- **If a fix doesn't change the crash (same address, same instruction, same
  backtrace), it didn't fix the actual bug** — even if it was a real and
  worthwhile fix for something else. Re-run the exact same reproduction
  after every fix rather than assuming a plausible-looking change worked;
  don't declare a crash fixed on code-review confidence alone.
- **Reach for a memory sanitizer (AddressSanitizer / `-fsanitize=address`, or
  Valgrind) as soon as a plain backtrace looks like a symptom rather than a
  cause** (a crash inside a destructor, in unrelated-looking code, or at a
  different spot each run). It reports the exact instruction and
  allocation/free site of the actual fault, which a post-hoc backtrace
  cannot reconstruct. An ASan build is usually just a separate build
  directory with `-fsanitize=address -fno-omit-frame-pointer` added to
  compile and link flags. For a qmake-based project, a build config that
  sets its defaults via unconditional `+=`/`-=` lines in the `.pro` file may
  silently re-override qmake command-line flag overrides (later assignments
  in the file, or file assignments processed after command-line ones, win);
  verify a debug/sanitizer build's *actual* effective flags in the generated
  Makefile before trusting it, rather than assuming the invocation worked.
- **Narrow the reproduction before diving into code.** Change one variable
  at a time between test runs ("does it still crash if I skip this one
  step?") to establish what's actually necessary to trigger the bug, rather
  than assuming every step in a reported sequence is essential.
- **Watch for self-referential container operations — but only on containers
  that can reallocate their elements.** `container.push_back(container
  .back())` (or any call passing a reference/pointer into a container as an
  argument to a method on that same container) is only a hazard for
  containers whose growth may move already-stored elements to a new backing
  buffer (e.g. `std::vector`/similar contiguous-storage containers). It is
  *not* an issue for containers that never relocate existing elements on
  insertion (e.g. `std::list`, `std::deque`'s existing elements, node-based
  maps/sets) — know which kind of container is involved before assuming a
  bug. Where it does apply, there are two valid fixes:
  - Copy the referenced value to a local before the call, so the argument no
    longer aliases storage the call might move. Simple and safe, but adds a
    copy on every call — avoid it in a hot path if the element type is large
    or expensive to copy.
  - Ensure the call cannot reallocate in the first place (e.g. `reserve()`
    enough capacity up front). This avoids the copy, but means every
    reference/pointer/iterator into the container — not just the one in this
    call — stays valid only as long as that reserved capacity isn't
    exceeded; a later change that adds one more element than planned for
    silently reintroduces the exact same class of bug. Prefer this only when
    the container's maximum size is genuinely known and enforced, and say so
    at the point capacity is reserved.
- **Object-ownership boundaries between manually-managed and
  garbage-collected/reference-counted systems are a common source of
  premature frees.** When a native object is exposed to a scripting engine,
  plugin system, or other embedded runtime, check that runtime's default
  ownership rule for a wrapped object (e.g. "no parent/owner set means the
  runtime takes ownership"). An externally-owned, long-lived object with no
  obvious owner marker can be silently deleted once the embedding runtime's
  GC (or refcounting) considers it unreachable — long after the code that
  exposed it has returned.
- **A single reported symptom can have multiple independent causes.** Don't
  stop at the first bug a sanitizer finds if the reproduction still crashes
  afterward — fix it, rebuild, and reproduce again, repeating until the
  actual reported symptom is gone. Keep unrelated fixes found along the way
  in their own commits rather than bundling them with the one that actually
  explains the report.

## General rules for development on this project

These were established while migrating the codebase from Qt5 to Qt6 and
fixing the regressions that surfaced, but apply broadly, not just to Qt work:

- **Prefer minimal, targeted changes.** Don't refactor, rename, or "clean up"
  code you're not there to touch as a side effect of an unrelated fix. If a
  change is inherently mechanical (a rename, an API replacement) and touches
  many files, keep each occurrence's diff as small as the mechanical change
  itself, not an opportunity for other improvements.
- **Comments describe current code only, tersely.** Don't explain what code
  used to do or why it was changed in an inline comment — that belongs in
  the commit message, which has room for it. The one exception: a short note
  is warranted when something looks non-standard specifically *because* it's
  preserving a behavior or working around a library limitation (a genuine
  "this looks wrong but isn't, because X" case) — keep even those short and
  focused on the current code, not a history lesson.
- **Commit messages carry the "why" and "what changed."** Explain the
  reasoning, the symptom being fixed, and any non-obvious consequences.
  Group commits by logical change, not by chronology — a mechanical rename
  across many files is one commit; a rewrite of one subsystem is one commit;
  don't bundle unrelated fixes together just because they happened in the
  same session. End commit messages with an attribution line crediting AI
  assistance where applicable.
- **Verify, don't assume — especially about library/framework behavior.**
  When a bug's cause hinges on how a library actually behaves (not just what
  its docs say), check: read the actual header/source if available, or write
  a small standalone program to test the specific behavior in question,
  rather than reasoning from memory alone. This caught real, non-obvious bugs
  during this migration (e.g. `QOpenGLWidget::makeCurrent()`'s FBO-rebinding
  side effect) that would have been easy to get wrong by assumption.
- **Rebuild after every change before considering it done.** A clean compile
  is cheap to verify and catches an enormous fraction of mistakes; there's no
  excuse to skip it even for a one-line change.
- **Be explicit about what you couldn't verify.** If interactive/GUI testing
  isn't possible in the current environment, say so plainly rather than
  implying a fix is confirmed. Static analysis and a clean build are not the
  same as confirming the actual reported symptom is gone.
- **Investigate unfamiliar state before acting on it, especially anything
  touching shared/external state (git branches, permissions, other running
  sessions).** Stop and ask rather than forcing through a permission error
  or an unexpected repository state.
- **When splitting work into commits, split by hunk/line when a single file's
  diff mixes unrelated concerns** — don't let commit granularity be
  constrained by which files happen to contain which changes.
