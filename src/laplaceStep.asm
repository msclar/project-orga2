%include "copiarFila.asm"

section .data
    format: db "num: %d", 10, 0
    mask4: DQ 4.0, 4.0, 4.0, 4.0

section .text
	global laplaceStep

laplaceStep:
	push rbp
	mov rbp, rsp
	push r12
	push r13
	push r14
	push r15
	
	%define res rdi
	%define phi rsi
	%define max_i r8
	%define max_j rcx
	%define s r10
	%define max_ij r11 ; por ahora no lo uso para nada
	%define i r12

	mov max_i, rdx
	mov s, max_j
	inc s ; s = max_j + 1 es la primer posicion de la matriz a analizar
	shl s, 3 ; porque cada posicion es de ocho bytes (64 bits)
	
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
	
	loop_i:
		%define sPhi r13
		
		mov sPhi, s ; se para en el medio de la estrella
		sub sPhi, max_j ; muevo una fila para arriba, al vecino #0
		vmovupd ymm1, [phi + sPhi] ; ymm1 es el acumulador
		
		add sPhi, max_j
		sub sPhi, 8
		vaddpd ymm1, [phi + sPhi]
		
		add sPhi, 8
		add sPhi, max_j
		vaddpd ymm1, [phi + sPhi]
		
		sub sPhi, max_j
		add sPhi, 8
		vaddpd ymm1, [phi + sPhi]
				
		vdivpd ymm1, [mask4]
		vmovupd [res + s], ymm1
				
		add s, 32 ; me muevo 4 doubles para la derecha

		dec i
		cmp i, 0
		jne loop_i

	mov r9, max_j
	call copiarbordes
	; fin choreo
	pop r15
	pop r14
	pop r13
	pop r12
	pop rbp
ret
