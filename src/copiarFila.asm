section .text
copiarfila: ;copia a la fila2 el contenido de la fila1
	%define input rdi
	%define fila1 r14
	%define fila2 r12
	%define max_j r9
	%define max_ij r11
	%define max_i r8
	%define s r10
	%define j r13		
	%define i r15
	mov j, max_j
	shr j, 2
	mov s, 0
	
	loop_fila:
		vmovupd ymm1, [fila1+s]
		vmovupd [fila2+s], ymm1
		add s, 32
		dec j
		cmp j, 0
		jne loop_fila
	
	mov s, max_j
	shl s, 3
	sub s, 32
	vmovupd ymm1, [fila1+s]
	vmovupd [fila2+s], ymm1
	ret

copiarbordes:
	mov s, max_j
	mov i, max_i
	sub i, 2
	loop_borde_izq:
		mov j, [input + s + 8]
		mov [input + s], j
		add s, max_j
		dec i
		cmp i, 0
		jne loop_borde_izq
	
	mov s, max_j
	add s, max_j
	sub s, 16
	mov i, max_i
	sub i, 2
	loop_borde_der:
		mov j, [input + s]
		mov [input + s + 8], j
		add s, max_j
		dec i
		cmp i, 0
		jne loop_borde_der
	
	%define fila1 r14
	%define fila2 r12
	
	mov fila2, input
	mov fila1, input
	add fila1, max_j
	shr max_j, 3
	call copiarfila
	
	mov fila2, input
	add fila2, max_ij
	shl max_j, 3
	sub fila2, max_j
	mov fila1, fila2
	sub fila1, max_j
	shr max_j, 3
	call copiarfila
	ret
