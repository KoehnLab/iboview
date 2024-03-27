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
    qtbase5-dev \
    qtscript5-dev \
    libqt5svg5-dev \
    libgl1-mesa-dev
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
    libqt5-qtbase-common-devel \
    libqt5-qtbase-devel \
    libqt5-qtsvg-devel \
    libqt5-qtscript-devel \
    Mesa-libGL-devel
```

### Compiling

Before starting, verify that `qmake --version` informs you that you are using Qt in version 5.x

```bash
mkdir build && cd build

# Select qmake or qmake-qt5 if qmake is not available on your system
QMAKE="qmake"
if [[ ! -x "$( command -v "$QMAKE" )" ]]; then
    QMAKE="qmake-qt5"
fi

$QMAKE ../main.pro
make -j $(nproc)
```

Afterwards, you'll have the `iboview` executable inside your `build` directory.
