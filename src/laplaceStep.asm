;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0

section .text
	global laplaceStep
    extern printf

laplaceStep:
	push rbp
	mov rbp, rsp
	push r14
	push r15
	
	%define res rdi
	%define phi rsi
	%define max_i rdx ; r8
	%define max_j rcx ; r9

		
	pop r15
	pop r14
	pop rbp
ret
