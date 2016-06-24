%include "copiarFila.asm"

section .data
    format: db "num: %d" , 10, 0
    fours: dq 4.0, 4.0, 4.0, 4.0
	eights: dq 8.0, 8.0, 8.0, 8.0

section .text
	global createA
    extern printf

createA:
	push rbp
	mov rbp, rsp
	push r12
	push r13
	push r14
	push r15

	%define A     rdi
	%define k     rsi
	%define max_i r8
	%define max_j rcx
	mov r12, A
	mov max_i, rdx
	
	%define indep_term xmm0
	%define delta_x    xmm1 
	%define delta_y    xmm2

	%define s r15
	mov s, max_j
	inc s
	shl s, 3 ; s = max_j + 1 (*8 bytes)

	%define max_ij r11
	mov rax, max_j
	mul max_i ; aqui max_i se destruye porque estaba en rdx, pero no lo uso mas
	mov max_ij, rax ; max_ij = max_i * max_j

	%define i r10
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

	add A, s ; ubico a A en la posicion s
	add k, s
	sub k, max_j ; lo dejo una fila arriba de s

	%define mask_indep_term ymm0 ; funciona porque indep_term esta en xmm0
	movddup indep_term, indep_term
	vinsertf128 mask_indep_term, mask_indep_term, indep_term, 1

	%define mask_delta_x_squared ymm1 ; funciona porque delta_x esta en xmm1	
	movddup delta_x, delta_x
	mulpd delta_x, delta_x
	vinsertf128 mask_delta_x_squared, mask_delta_x_squared, delta_x, 1

	%define mask_delta_y_squared ymm2 ; funciona porque delta_y esta en xmm2
	movddup delta_y, delta_y
	mulpd delta_y, delta_y
	vinsertf128 mask_delta_y_squared, mask_delta_y_squared, delta_y, 1
	
	%define mask_4_times_delta_x_squared ymm7
	vmovupd mask_4_times_delta_x_squared, [fours]
	vmulpd mask_4_times_delta_x_squared, mask_delta_x_squared

	%define mask_4_times_delta_y_squared ymm8
	vmovupd mask_4_times_delta_y_squared, [fours]
	vmulpd mask_4_times_delta_y_squared, mask_delta_y_squared
	
	loop_i:
		;~ haciendo un poco de aritmetica llegamos a lo siguiente:
		;~ A[indice(0, s, max_i*max_j)] = (- k[indice(i+1, j, max_j)] + k[indice(i-1, j, max_j)] + 4 * k[indice(i-1, j, max_j)]) / (4 * delta_x * delta_x);
		;~ A[indice(2, s, max_i*max_j)] = (k[indice(i+1, j, max_j)] - k[indice(i-1, j, max_j)] + 4 * k[indice(i+1, j, max_j)]) / (4 * delta_x * delta_x);

		vmovupd ymm3, [k]             ; k[indice(i-1, j, max_j)]
		vmovupd ymm4, [k + 2 * max_j] ; k[indice(i+1, j, max_j)]
		;~ vmovupd ymm5, [k + max_j]     ; k[indice(i, j, max_j)]

		; procesamos el vecino #0 de A		
		vmovupd ymm6, ymm3
		vmulpd ymm6, [fours] ; 4 * k[indice(i-1, j, max_j)]
		vsubpd ymm6, ymm4     ; - k[indice(i+1, j, max_j)] + 4 * k[indice(i-1, j, max_j)]
		vaddpd ymm6, ymm3     ; - k[indice(i+1, j, max_j)] + k[indice(i-1, j, max_j)] + 4 * k[indice(i-1, j, max_j)]
		vdivpd ymm6, mask_4_times_delta_x_squared
		vmovupd [A], ymm6
		
		; procesamos el vecino #2 de A		
		vmovupd ymm6, ymm4
		vmulpd ymm6, [fours] ; 4 * k[indice(i+1, j, max_j)]
		vsubpd ymm6, ymm3     ; - k[indice(i-1, j, max_j)] + 4 * k[indice(i+1, j, max_j)]
		vaddpd ymm6, ymm4     ; - k[indice(i-1, j, max_j)] + k[indice(i+1, j, max_j)] + 4 * k[indice(i+1, j, max_j)]
		vdivpd ymm6, mask_4_times_delta_x_squared
		vmovupd [A + 2 * max_ij], ymm6
		
		
		;~ A[indice(1, s, max_i*max_j)] = (- k[indice(i, j+1, max_j)] + k[indice(i, j-1, max_j)] + 4 * k[indice(i, j-1, max_j)]) / (4 * delta_y * delta_y)
		;~ A[indice(3, s, max_i*max_j)] = (  k[indice(i, j+1, max_j)] - k[indice(i, j-1, max_j)] + 4 * k[indice(i, j+1, max_j)]) / (4 * delta_y * delta_y)

		
		vmovupd ymm3, [k + max_j - 8] ; k[indice(i, j-1, max_j)]
		vmovupd ymm4, [k + max_j + 8] ; k[indice(i, j+1, max_j)]
		
		; procesamos el vecino #1 de A		
		vmovupd ymm6, ymm3
		vmulpd ymm6, [fours] ; 4 * k[indice(i, j-1, max_j)]
		vsubpd ymm6, ymm4    ; - k[indice(i, j+1, max_j)] + 4 * k[indice(i, j-1, max_j)]
		vaddpd ymm6, ymm3    ; - k[indice(i, j+1, max_j)] + k[indice(i, j-1, max_j)] + 4 * k[indice(i-1, j, max_j)]
		vdivpd ymm6, mask_4_times_delta_y_squared
		vmovupd [A + max_ij], ymm6
		
		; procesamos el vecino #3 de A		
		vmovupd ymm6, ymm4
		vmulpd ymm6, [fours] ; 4 * k[indice(i, j+1, max_j)]
		vsubpd ymm6, ymm3    ; - k[indice(i, j-1, max_j)] + 4 * k[indice(i, j+1, max_j)]
		vaddpd ymm6, ymm4    ; - k[indice(i, j-1, max_j)] + k[indice(i, j+1, max_j)] + 4 * k[indice(i, j+1, max_j)]
		vdivpd ymm6, mask_4_times_delta_y_squared
		vmovupd [A + r14], ymm6 ; r14 = 3 * max_ij
		
		;~ A[indice(4, s, max_i*max_j)] = - 8 * k[s] / (4 * delta_x * delta_x) - 8 * k[s] / (4 * delta_y * delta_y) - w_b * C_b * rho_b - rho * C_rho / delta_t;
		
		vmovupd ymm3, [k + max_j] ; k[s]
		vmulpd ymm3, [eights]
		vmovapd ymm4, ymm3
		
		vdivpd ymm3, mask_4_times_delta_x_squared ; - 8 * k[s] / (4 * delta_x * delta_x)
		vdivpd ymm4, mask_4_times_delta_y_squared ; - 8 * k[s] / (4 * delta_y * delta_y)
		
		vmovapd ymm5, mask_indep_term
		vsubpd ymm5, ymm3
		vsubpd ymm5, ymm4
		vmovupd [A + 4 * max_ij], ymm5
		
		add A, 32
		add k, 32
		
		dec i
		jnz loop_i
	
	mov A, r12
	mov rsi, 5
	
	loop_bordes:
		mov r9, max_j
		call borrarbordes
		add A, max_ij
		dec rsi
		jnz loop_bordes
	
	pop r15
	pop r14
	pop r13
	pop r12
	pop rbp
ret
