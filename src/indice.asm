global indice
section .text

indice2:
	push rbp
	mov rbp, rsp

	%define i rdi
	%define j rsi
	%define max_j rdx
	
	mov rax, i
	mul max_j
	add rax, j
	
	pop rbp
ret


; nasm -f elf64 funcion.asm -o funcion.o
; gcc -o ejec programa.c funcion.o
