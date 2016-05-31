CC = gcc
ASM = nasm
DBG = gdb
CFLAGS = -ggdb -Wall -std=c99 -pedantic -m64 -lm
ASMFLAGS = -felf64 -g -F dwarf

BIN = jacobi-optimizado
BIN_DIR = ../bin

OBJS = jacobi-optimizado.o indice.o \
			jacobiStep.o

.PHONY: all clean

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm
	cp $(BIN) $(BIN_DIR)/$(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.asm
	$(ASM) $(ASMFLAGS) -o $@ $<

all: $(BIN)

clean:
	rm -f $(OBJS)
	rm -f ./$(BIN)
	rm -f $(BIN_DIR)/$(BIN)
