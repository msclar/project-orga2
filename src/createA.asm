section .data
    format: db "num: %d" , 10, 0
    zeros: dq 0.0, 0.0

section .text
	global createA
    extern printf

createA:
	push rbp
	mov rbp, rsp

	%define A     rdi
	%define k     rsi
	%define max_i rdx
	%define max_j rcx
	
	%define indep_term xmm0
	%define delta_x    xmm1 
	%define delta_y    xmm2

	pop rbp
ret
