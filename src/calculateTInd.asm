;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0
    mask4: DQ 4.0, 4.0, 4.0, 4.0
	twos: DQ 2.0, 2.0

section .text
	global calculateTInd
    extern printf

calculateTInd:
	push rbp
	mov rbp, rsp
	push r12
	push r13
	
	%define phi rdi
	%define sigma rsi
	%define TInd r9
	%define max_i rcx
	%define max_j r8
	%define delta_x xmm0
	%define delta_y xmm1
	%define resto xmm2
	%define s r10
	%define max_ij r11
	%define i r12
	%define j r13
	
	%define mask_two_delta_y ymm1 ; funciona porque delta_y esta en xmm1
	mulsd delta_y, [twos]
	movddup delta_y, delta_y
	vinsertf128 mask_two_delta_y, mask_two_delta_y, delta_y, 1
	
	%define mask_two_delta_x ymm0 ; funciona porque delta_y esta en xmm1
	mulsd delta_x, [twos]
	movddup delta_x, delta_x
	vinsertf128 mask_two_delta_x, mask_two_delta_x, delta_x, 1
	
	%define mask_resto ymm2 ; funciona porque delta_y esta en xmm1
	movddup resto, resto
	vinsertf128 mask_resto, mask_resto, resto, 1
	
	mov TInd, rdx

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
	
	add phi, 8
	loop_i:
		vmovupd ymm3, [phi + 2*max_j]
		vsubpd ymm3, [phi]
		vdivpd ymm3, mask_two_delta_x
		
		vmovupd ymm4, [phi + max_j + 8]
		vsubpd ymm4, [phi + max_j - 8]
		vdivpd ymm4, mask_two_delta_y
		
		vmulpd ymm3, ymm3
		vmulpd ymm4, ymm4
		vaddpd ymm3, ymm4
		vmulpd ymm3, [sigma + s]
		vaddpd ymm3, mask_resto
		
		vmovupd [TInd + s], ymm3
				
		add phi, 32
		add s, 32 ; me muevo 4 doubles para la derecha

		dec i
		cmp i, 0
		jne loop_i
		
	pop r13
	pop r12
	pop rbp
ret
