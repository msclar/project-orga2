%include "copiarFila.asm"

section .data
    format: db "num: %d", 10, 0

section .text
	global calculateVectorError

calculateVectorError:
	push rbp
	mov rbp, rsp
	push r12
	push r13
	push r14
	push r15
	
	%define A     rdi
	%define Tn    rsi
	%define B     rdx
	%define res   rcx
	%define max_i r8
	%define max_j r9
	mov r15, res
	
	%define s r10
	mov s, max_j
	inc s
	shl s, 3 ; s = max_j + 1 (*8 bytes)
	
	add A, s ; ubico a A en la posicion s
	add B, s
	add Tn, s
	add res, s
	
	mov r13, B

	%define max_ij r11
	mov rax, max_j
	mul max_i
	mov max_ij, rax ; max_ij = max_i * max_j

	mov B, r13

	%define i r12
	mov i, max_ij
	sub i, max_j
	sub i, max_j
	inc i
	shr i, 2 ; i es la cantidad de casillas que importan (ceil ((max_ij - 2 * max_j - 2) / 4))
			 ; i = floor ((max_ij - 2 * max_j + 1) / 4)

	shl max_ij, 3
	shl max_j, 3
	
	mov r14, max_ij
	add r14, max_ij
	add r14, max_ij ; r14 = 3 * max_ij

	sub Tn, max_j ; lo coloco una fila arriba de lo que tendria que estar para las cuentas

	loop_i:
		vxorpd ymm3, ymm3
		vsubpd ymm3, [B] ; ymm3 = -B
	
		vmovupd ymm2, [Tn] ; cargo el vecino de arriba y los 3 elementos que siguen en Tn
		vmulpd ymm2, [A]
		vaddpd ymm3, ymm2 ; ymm3 es el acumulador
		
		; cargo el vecino de la izquierda 
		vmovupd ymm2, [Tn + max_j - 8]
		vmulpd ymm2, [A + max_ij]
		vaddpd ymm3, ymm2
		
		; cargo el vecino de abajo
		vmovupd ymm2, [Tn + 2 * max_j]
		vmulpd ymm2, [A + 2 * max_ij]
		vaddpd ymm3, ymm2
		
		; cargo el vecino de la derecha
		vmovupd ymm2, [Tn + max_j + 8]
		vmulpd ymm2, [A + r14] ; r14 = 3 * max_ij
		vaddpd ymm3, ymm2
		
		; cargo el centro
		vmovupd ymm2, [Tn + max_j]
		vmulpd ymm2, [A + 4 * max_ij]
		vaddpd ymm3, ymm2
		
		; cargo los resultados
		vmovupd [res], ymm3
		
		add A, 32
		add B, 32
		add Tn, 32
		add res, 32
		
		dec i
		jnz loop_i
	
	mov rdi, r15
	call borrarbordes

	pop r15
	pop r14
	pop r13
	pop r12
	pop rbp
ret

