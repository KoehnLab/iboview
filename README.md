# IboView

Source code of the IboView program (http://iboview.org/). Unfortunately, the original program is no longer maintained and the source code as provided
on their website fails to compile with modern compilers that actually require C++ standard-conforming code (which IboView in parts isn't). Therefore,
this repo contains our own "fork" of the project in order to incorporate necessary patches ourselves (and potentially also some quality-of-life
improvements).

## Building

### Installing dependencies

#### Ubuntu

```bash
sudo apt-get install \
    build-essential \
    libboost-all-dev \
    qt6-base-dev \
    qt6-declarative-dev \
    libqt6svg6-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev
```

#### OpenSuse

```bash
# This is only to figure out which version of Boost to install as OpenSuse doesn't
# seem to have an unversioned default Boost package but has different versions in the repos
BOOST_PKGS="$( zypper search libboost_atomic*-devel | grep "libboost" | cut -d "|" -f 2 )"
BOOST_VERSION=$( echo "$BOOST_PKGS" | sort | tail -n 1 | tr -cd '[0-9_\n]' | sed 's/^_\+//' )
echo "Using Boost version $BOOST_VERSION"
BOOST_PKGS="$( zypper search libboost_*-devel | grep "$BOOST_VERSION" | cut -d "|" -f 2 | tr '\n' ' ' )"

sudo zypper install \
    gcc \
    make \
    glu-devel \
    $BOOST_PKGS \
    qt6-base-common-devel \
    qt6-base-devel \
    qt6-svg-devel \
    qt6-declarative-devel \
    Mesa-libGL-devel
```

#### Arch-Based

```bash
yay  -S qt6-base \
        openblas \
        intel-mkl \
        qt6-declarative \
        glu
```

### Compiling

Before starting, verify that `qmake6 --version` (or `qmake --version`, if
`qmake6` isn't available under that name on your system) informs you that
you are using Qt in version 6.x. Prefer `qmake6` when both exist: on some
systems (e.g. Ubuntu 22.04) a generic, unversioned `qmake` exists only as a
`qtchooser` wrapper with no Qt6 profile registered, and fails outright
rather than falling back to Qt6.

```bash
mkdir build && cd build

# Prefer qmake6, falling back to plain qmake only if qmake6 isn't installed
# under that name. Don't do this the other way around: an unversioned
# "qmake" found on PATH may be a non-functional qtchooser wrapper that has
# no Qt6 profile configured, even though it exists and is executable.
QMAKE="qmake6"
if [[ ! -x "$( command -v "$QMAKE" )" ]]; then
    QMAKE="qmake"
fi

$QMAKE ../main.pro
make -j $(nproc)
```

Afterwards, you'll have the `iboview` executable inside your `build` directory.
