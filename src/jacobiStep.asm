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
	mov rbp, rsp
	push r14
	push r15
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
	
	mov i, max_i ; esto seguro esta mal
	sub i, 2
	mov s, max_j
	add s, 1
	shl s, 3
	
	mov max_ij, max_i
	mov rax, max_j
	mul max_ij
	mov max_ij, rax
	shl max_ij, 3
	
	loop_i:
		mov j, max_j
		shl max_j, 3
		add j, 1
		shr j, 2
		loop_j:
			;vmovupd: leer 4 doubles contiguos de memoria
			%define sTn r14
			mov sTn, s
			sub sTn, max_j
			vxorpd ymm3, ymm3
			vmovupd ymm1, [Tn + sTn]
			vmovupd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add sTn, max_j
			sub sTn, 8
			add A, max_ij
			vmovupd ymm1, [Tn + sTn]
			vmovupd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add sTn, 8
			add sTn, max_j
			add A, max_ij
			vmovupd ymm1, [Tn + sTn] 
			vmovupd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			sub sTn, max_j
			add sTn, 8
			add A, max_ij
			vmovupd ymm1, [Tn + sTn]
			vmovupd ymm2, [A + s]
			vmulpd ymm2, ymm1
			vaddpd ymm3, ymm2
			
			add A, max_ij
			vmovupd ymm1, [B + s]
			vmovupd ymm2, [A + s]
			vsubpd ymm1, ymm3
			vdivpd ymm1, ymm2
			vmovupd [Tn_sig + s], ymm1 
			
			mov rax, max_ij
			shl rax, 2
			sub A, rax
			
			add s, 32
			dec j
			cmp j, 0
			jne loop_j
		
		shr max_j, 3
		dec i
		cmp i, 0
		jne loop_i
		
	pop r15
	pop r14
	pop rbp
ret
