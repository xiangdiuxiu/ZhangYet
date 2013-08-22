all:src/main.cpp src/main.h src/readBinData.cpp src/readBinData.h src/matrix.cpp src/matrix.h src/model.cpp src/model.h src/Individual.cpp src/Individual.h src/Locus.cpp src/Locus.h
	g++ src/main.cpp src/main.h src/readBinData.cpp src/readBinData.h src/matrix.cpp src/matrix.h src/model.cpp src/model.h src/Individual.cpp src/Individual.h src/Locus.cpp src/Locus.h -o mplink

install:
	-mkdir -p $(RPM_INSTALL_ROOT)/usr/local/bin/
	install -m 755 mplink $(RPM_INSTALL_ROOT)/usr/local/bin/mplink