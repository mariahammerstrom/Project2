TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo
LIBS += -lunittest++

SOURCES += main.cpp \
    lib.cpp \
    tests.cpp

HEADERS += \
    lib.h

