TEMPLATE = app
TARGET = boost_special_test
DEPENDPATH += .

CONFIG += c++11
QT =
CONFIG -= debug
CONFIG += release

QMAKE_CXXFLAGS += -Wno-deprecated-copy -Wno-class-memaccess -Wno-implicit-fallthrough

QMAKE_CXXFLAGS_RELEASE -= -O2 -mtune=generic
QMAKE_CXXFLAGS_RELEASE -= -fstack-protector -fstack-protector-strong --param=ssp-buffer-size=4
QMAKE_CXXFLAGS_RELEASE += -Ofast -ffast-math -march=native -DNDEBUG -DINCLUDE_OPTIONALS

SOURCES += main.cpp
