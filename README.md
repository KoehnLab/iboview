# IboView

Source code of the IboView program (http://iboview.org/). Unfortunately, the original program is no longer maintained and the source code as provided
on their website fails to compile with modern compilers that actually require C++ standard-conforming code (which IboView in parts isn't). Therefore,
this repo contains our own "fork" of the project in order to incorporate necessary patches ourselves (and potentially also some quality-of-life
improvements).

## Building

### Installing dependencies

#### Ubuntu

```bash
sudo apt install \
    libboost-all-dev \
    qtbase5-dev \
    qtscript5-dev \
    libqt5svg5-dev \
    libgl1-mesa-dev
```

### Compiling

Before starting, verify that `qmake --version` informs you that you are using Qt in version 5.x

```bash
mkdir build && cd build
qmake ../main.pro
make -j $(nproc)
```
