;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0
    zeros: dq 0.0, 0.0

section .text
	global fillWithZeros
    extern printf

fillWithZeros:
	push rbp
	mov rbp, rsp
	
	%define matrix          rdi
	%define max_ij          rsi

	%define all_zeros ymm0
	vmovupd xmm0, [zeros]
	vinsertf128 all_zeros, all_zeros, xmm0, 1

	%define i r8
	mov i, max_ij
	sub i, 4 ; i = max_ij - 4
	
	; proceso los ultimos 4 elementos
	vmovupd [matrix + 8*i], all_zeros
	
	; redondeo para abajo a un multiplo de 4 (idem updateB.asm)
	dec i
	shr i, 2
	shl i, 5

	add i, 32 ; para contrarrestar lo que voy a restar despues
	loop_i:
		sub i, 32 ; saco los 4 elementos ya analizados (*8 bytes)
		
		vmovupd [matrix + i], all_zeros
		jnz loop_i 
	
	pop rbp
ret
