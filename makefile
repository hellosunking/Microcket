all: bin/krmdup bin/krmdup.pipe bin/sam2pairs

GXXFlag=-std=c++11 -O3

bin/krmdup: src/preprocess/krmdup.cpp
	@echo Build Krmdup
	@cd src/preprocess/; g++ krmdup.cpp -fopenmp -o ../../bin/krmdup $(GXXFlag) -lz ; cd ../../

bin/krmdup.pipe: src/preprocess/krmdup.pipe.cpp
	@echo Build Krmdup
	@cd src/preprocess/; g++ krmdup.pipe.cpp -fopenmp -o ../../bin/krmdup.pipe $(GXXFlag) -lz ; cd ../../

bin/sam2pairs: src/sam2pairs/pairutil.h src/sam2pairs/flash2pairs.h src/sam2pairs/unc2pairs.h src/sam2pairs/sam2pairs.cpp
	@echo Build sam2pairs
	@cd src/sam2pairs/; g++ sam2pairs.cpp -fopenmp -o ../../bin/sam2pairs $(GXXFlag) -lz ; cd ../../

clean:
	@rm -f bin/ktrim.mod bin/krmdup bin/cut.fq.tail bin/sam2pairs

