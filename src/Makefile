CC = g++
CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -fopenmp -march=native -mavx2
BUILD_DIR = build

TARGET = ${BUILD_DIR}/Main
OBJS = ${BUILD_DIR}/Main.o \
       ${BUILD_DIR}/BM.o \
       ${BUILD_DIR}/BP.o \
       ${BUILD_DIR}/KMP.o \
       ${BUILD_DIR}/HybridPicker.o \
       ${BUILD_DIR}/FastaReader.o \
       ${BUILD_DIR}/Benchmark.o

# --- Linking step ---
${TARGET}: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -o ${TARGET}

# --- Compilation steps ---
${BUILD_DIR}/Main.o: main.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/BM.o: imp/BMh.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/BP.o: imp/BPh.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/KMP.o: imp/KMPh.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/HybridPicker.o: imp/HybridPicker.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/FastaReader.o: imp/FastaReader.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}/Benchmark.o: imp/Benchmark.cpp | ${BUILD_DIR}
	${CC} ${CFLAGS} -c $< -o $@

${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}

clean:
	rm -rf ${BUILD_DIR}

.PHONY: clean
