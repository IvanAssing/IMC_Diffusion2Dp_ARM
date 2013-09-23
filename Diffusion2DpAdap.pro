TEMPLATE = app
CONFIG += console
#CONFIG -= app_bundle
#CONFIG -= qt

QT       += core gui opengl
QT += widgets

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS *= -fopenmp

LIBS += -lquadmath

SOURCES += main.cpp \
    node2d.cpp \
    imc_dfm.cpp \
    functor2d.cpp \
    boundary2d.cpp \
    data.cpp \
    gaussseidel.cpp \
    graphics.cpp \
    diffusion2dpar.cpp \
    element2d.cpp \
    solverlu.cpp

HEADERS += \
    node2d.h \
    imc_dfm.h \
    functor2d.h \
    boundary2d.h \
    data.h \
    gaussseidel.h \
    graphics.h \
    diffusion2dpar.h \
    element2d.h \
    solverlu.h

FORMS +=

