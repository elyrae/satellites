PROJECT_ROOT_PATH = $${PWD}/

win32: OS_SUFFIX = win32
linux-g++: OS_SUFFIX = linux

QMAKE_LIBS+= -lgomp -lpthread

CONFIG(debug, debug|release) {
    BUILD_FLAG = debug
    LIB_SUFFIX = d
} else {
    BUILD_FLAG = release
    QMAKE_CXXFLAGS += -fopenmp -march=native -O2 -finline-functions -ftree-vectorize -mavx -mavx2
}

# Sources
LIBS_PATH   = $${PROJECT_ROOT_PATH}/lib.$${OS_SUFFIX}/
INC_PATH    = $${PROJECT_ROOT_PATH}/include/
IMPORT_PATH = $${PROJECT_ROOT_PATH}/import/
BIN_PATH    = $${PROJECT_ROOT_PATH}/bin/$${BUILD_FLAG}/
SRC_PATH    = $${PROJECT_ROOT_PATH}/src/

# Builds paths
BUILD_PATH = $${PROJECT_ROOT_PATH}/build/$${BUILD_FLAG}/$${TARGET}/
RCC_DIR     = $${BUILD_PATH}/rcc/
UI_DIR      = $${BUILD_PATH}/ui/
MOC_DIR     = $${BUILD_PATH}/moc/
OBJECTS_DIR = $${BUILD_PATH}/obj/

LIBS += -L$${LIBS_PATH}/
# INCLUDEPATH += $${INC_PATH}/
# INCLUDEPATH += $${IMPORT_PATH}/
INCLUDEPATH = $${SRC_PATH}/include/

CONFIG += c++14

QMAKE_RPATHDIR += /home/elyrae/Qt/5.7/gcc_64/lib
