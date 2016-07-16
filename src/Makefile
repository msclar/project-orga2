CC = gcc
ASM = nasm
DBG = gdb
CFLAGS = -ggdb -DASM -DMEASURE_TIME -Wall -O2 -pedantic -m64 -lm
ASMFLAGS = -felf64 -g -F dwarf

BIN = jacobi-optimizado
BIN_DIR = 

OBJS = jacobi-optimizado.o fillWithConstant.o calculateVectorError.o \
       jacobiStep.o laplaceStep.o updateB.o createA.o \
       calculateTInd.o vectorNorm.o distanceBetweenVectors.o

.PHONY: all clean

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

%.o: %.c
	mkdir -p output && $(CC) $(CFLAGS) -c -o $@ $<

%.o: %.asm
	$(ASM) $(ASMFLAGS) -o $@ $<

all: $(BIN)

clean:
	rm -f $(OBJS)
	rm -f ./$(BIN)
	rm -f $(BIN_DIR)/$(BIN)
