;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0

section .text
	global jacobiStep
    extern printf

jacobiStep:
	push rbp
	push r14
	push r15
	mov rbp, rsp

	%define Tn_sig rdi
	%define Tn rsi
	%define B rdx
	%define A rcx
	%define max_i r8
	%define max_j r9
	
	%define s r10
	%define max_ij r11
	%define i r12
	%define j r13
	
	mov i, max_i ; esto seguro esta mal
	sub i, 2
	mov s, max_j
	add s, 1
	shl s, 4
	
	mov max_ij, max_i
	mov rax, max_j
	mul max_ij
	mov max_ij, rax
	shl max_ij, 4
	
	loop_i:
		mov j, max_j
		shl max_j, 4
		sub j, 2
		shr j, 2
		
		loop_j:
			;vmovapd: leer 4 doubles contiguos de memoria
			%define sTn r14
			mov sTn, s
			sub sTn, max_j
			vmovapd ymm1, [Tn + sTn]
			vmovapd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add sTn, max_j
			sub sTn, 16
			add A, max_ij
			vmovapd ymm1, [Tn + sTn]
			vmovapd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add sTn, 16
			add sTn, max_j
			add A, max_ij
			vmovapd ymm1, [Tn + sTn] 
			vmovapd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			sub sTn, max_j
			add sTn, 16
			add A, max_ij
			vmovapd ymm1, [Tn + sTn]
			vmovapd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add A, max_ij
			vmovapd ymm1, [B + s]
			vmovapd ymm2, [A + s]
			vsubpd ymm1, ymm3
			vdivpd ymm1, ymm2
			vmovapd ymm1, [Tn_sig + s]
			
			mov rax, max_ij
			shl rax, 2
			sub A, rax
			
			add s, 64
			dec j
			cmp j, 0
			jne loop_j
		
		shr max_j, 4
		add s, 32
		dec i
		cmp i, 0
		jne loop_i
		
	pop r15
	pop r14
	pop rbp
ret
