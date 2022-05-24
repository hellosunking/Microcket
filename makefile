all: bin/ktrim.mod bin/krmdup bin/cut.fq.tail bin/sam2pairs

GXXFlag=-march=native -std=c++11 -O2

bin/ktrim.mod: src/preprocess/ktrim.cpp src/preprocess/common.h src/preprocess/util.h src/preprocess/param_handler.h src/preprocess/pe_handler.h
	@echo Build Ktrim
	@cd src/preprocess/; g++ -fopenmp -lz $(GXXFlag) -o ../../bin/ktrim.mod ktrim.cpp; cd ../../

bin/krmdup: src/preprocess/krmdup.cpp
	@echo Build Krmdup
	@cd src/preprocess/; g++ -fopenmp $(GXXFlag) -o ../../bin/krmdup krmdup.cpp; cd ../../

bin/cut.fq.tail: src/preprocess/cut.fq.tail.cpp
	@echo Build cut.fq.tail
	@cd src/preprocess/; g++ -fopenmp $(GXXFlag) -o ../../bin/cut.fq.tail cut.fq.tail.cpp; cd ../../

bin/sam2pairs: src/sam2pairs/pairutil.h src/sam2pairs/flash2pairs.h src/sam2pairs/unc2pairs.h src/sam2pairs/sam2pairs.cpp
	@echo Build sam2pairs
	@cd src/sam2pairs/; g++ -fopenmp $(GXXFlag) -o ../../bin/sam2pairs sam2pairs.cpp; cd ../../

clean:
	@rm -f bin/ktrim.mod bin/krmdup.mod bin/sam2pair

