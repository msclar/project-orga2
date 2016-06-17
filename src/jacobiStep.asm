;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0

section .text
	global jacobiStep
    extern printf

	%define Tn_sig rdi
	%define Tn rsi
	%define B r15
	%define A rcx
	%define max_i r8
	%define max_j r9
	mov r15, rdx
	%define s r10
	%define max_ij r11
	%define i r12
	%define j r13
	
copiarfila: ;copia a la fila2 el contenido de la fila1
	%define fila1 r14
	%define fila2 r12
	
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
		

jacobiStep:
	push rbp
	mov rbp, rsp
	push r12
	push r13
	push r14
	push r15
	
	mov B, rdx
	
	mov s, max_j
	inc s ; s = max_j + 1 es la primer posicion de la matriz a analizar
	shl s, 3 ; porque cada posicion es de ocho bytes (64 bits)
	add A, s
	
	mov rax, max_j
	mul max_i
	mov max_ij, rax ; max_ij = max_i * max_j
	
	mov i, max_ij
	sub i, max_j
	sub i, max_j
	inc i
	shr i, 2 ; i es la cantidad de casillas que importan (ceil ((max_ij - 2 * max_j - 2) / 4))
			 ; i = floor ((max_ij - 2 * max_j + 1) / 4)
	
	shl max_ij, 3
	shl max_j, 3
	
	lea Tn, [Tn + 8]
	
	mov r14, max_ij
	add r14, max_ij
	add r14, max_ij
	
	loop_i:
		vmovupd ymm2, [Tn] ; cargo el vecino #0 y los 3 elementos que siguen en Tn
		vmulpd ymm2, [A]
		vmovapd ymm3, ymm2 ; ymm3 es el acumulador
		
		vmovupd ymm2, [Tn + max_j - 8]
		vmulpd ymm2, [A + max_ij]
		vaddpd ymm3, ymm2
		
		vmovupd ymm2, [Tn + 2 * max_j]
		vmulpd ymm2, [A + 2 * max_ij]
		vaddpd ymm3, ymm2
		
		vmovupd ymm2, [Tn + max_j + 8]
		vmulpd ymm2, [A + r14]
		vaddpd ymm3, ymm2
		
		; Tn_sig[s] = (B[s] - sum) / A[indice(4, s, max_i * max_j)];
		vmovupd ymm1, [B + s]
		vsubpd ymm1, ymm3
		vdivpd ymm1, [A + 4 * max_ij]
		vmovupd [Tn_sig + s], ymm1 
		
		add A, 32
		add Tn, 32
		add s, 32 ; me muevo 4 doubles para la derecha

		dec i
		cmp i, 0
		jne loop_i
	
	;~ mov s, max_j
	;~ mov i, max_i
	;~ sub i, 2
	;~ loop_borde_izq:
		;~ mov j, [Tn_sig + s + 8]
		;~ mov [Tn_sig + s], j
		;~ add s, max_j
		;~ dec i
		;~ cmp i, 0
		;~ jne loop_borde_izq
	;~ 
	;~ mov s, max_j
	;~ add s, max_j
	;~ sub s, 16
	;~ mov i, max_i
	;~ sub i, 2
	;~ loop_borde_der:
		;~ mov j, [Tn_sig + s]
		;~ mov [Tn_sig + s + 8], j
		;~ add s, max_j
		;~ dec i
		;~ cmp i, 0
		;~ jne loop_borde_der
	;~ 
	
	;~ mov fila2, Tn_sig
	;~ mov fila1, Tn_sig
	;~ add fila1, max_j
	;~ shr max_j, 3
	;~ call copiarfila
	;~ 
	;~ mov fila2, Tn_sig
	;~ add fila2, max_ij
	;~ shl max_j, 3
	;~ sub fila2, max_j
	;~ mov fila1, fila2
	;~ sub fila1, max_j
	;~ shr max_j, 3
	;~ call copiarfila
	
	pop r15
	pop r14
	pop r13
	pop r12
	pop rbp
ret
