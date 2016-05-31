;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila

global jacobiStep
section .text

jacobiStep:
	push rbp
	mov rbp, rsp

	%define Tn_sig rdi
	%define Tn rsi
	%define B rdx
	%define A rcx
	%define max_i r8d
	%define max_j r9d
	
	%define s r10
	%define i r12d
	%define j r13d
	
	mov i, max_i ; esto seguro esta mal
	sub i, 2
	
	
	loop_i:
		mov j, max_j
		sub j, 2
		shr j, 2
		loop_j:
			;vmovapd: leer 4 doubles contiguos de memoria
			vmovapd ymm0, [Tn + s + 16] ;ymm0 = Tn[s+1] 
			
			
			dec j
			cmp j, 0
			jne loop_j
		dec i
		cmp i, 0
		jne loop_i
		
	
	pop rbp
ret
