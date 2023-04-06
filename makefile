all: bin/ktrim.mod bin/krmdup bin/krmdup.append bin/krmdup.pipe bin/cut.fq.tail bin/sam2pairs

GXXFlag=-std=c++11 -O2

bin/ktrim.mod: src/preprocess/ktrim.mod.cpp src/preprocess/common.h src/preprocess/util.h src/preprocess/param_handler.h src/preprocess/pe_handler.h
	@echo Build Ktrim
	@cd src/preprocess/; g++ ktrim.mod.cpp -o ../../bin/ktrim.mod $(GXXFlag) -fopenmp -lz ; cd ../../

bin/krmdup: src/preprocess/krmdup.cpp
	@echo Build Krmdup
	@cd src/preprocess/; g++ krmdup.cpp -o ../../bin/krmdup -fopenmp $(GXXFlag) ; cd ../../

bin/krmdup.append: src/preprocess/krmdup.append.cpp
	@echo Build Krmdup.append
	@cd src/preprocess/; g++ krmdup.append.cpp -o ../../bin/krmdup.append -fopenmp $(GXXFlag) ; cd ../../

bin/krmdup.pipe: src/preprocess/krmdup.pipe.cpp
	@echo Build Krmdup.pipe
	@cd src/preprocess/; g++ krmdup.pipe.cpp -o ../../bin/krmdup.pipe -fopenmp $(GXXFlag) ; cd ../../

bin/cut.fq.tail: src/preprocess/cut.fq.tail.cpp
	@echo Build cut.fq.tail
	@cd src/preprocess/; g++ cut.fq.tail.cpp -o ../../bin/cut.fq.tail -fopenmp $(GXXFlag) ; cd ../../

bin/sam2pairs: src/sam2pairs/pairutil.h src/sam2pairs/flash2pairs.h src/sam2pairs/unc2pairs.h src/sam2pairs/sam2pairs.cpp
	@echo Build sam2pairs
	@cd src/sam2pairs/; g++ sam2pairs.cpp -o ../../bin/sam2pairs -fopenmp $(GXXFlag) ; cd ../../

clean:
	@rm -f bin/ktrim.mod bin/krmdup bin/cut.fq.tail bin/sam2pairs

