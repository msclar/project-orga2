;Una funcion para cumplir con la Convencion C debe:
;Preservar RBX, R12, R13, R14 y R15
;Retornar el resultado en RAX o XMM0
;No romper la pila
section .data
    format: db "num: %d" , 10, 0
    zeros: dq 0.0, 0.0

section .text
	global updateB
    extern printf

updateB:
	push rbp
	mov rbp, rsp
	
	%define B               rdi
	%define Tn              rsi
	%define TIndAct         rdx
	%define rho_times_C_rho xmm0
	%define delta_t         xmm1
	%define max_ij          rcx

	vmovupd xmm2, [zeros] ; xmm2 = 0
	vsubpd xmm2, rho_times_C_rho
	vmovapd rho_times_C_rho, xmm2 ; rho_times_C_rho = - rho * C_rho
	
	%define mask_rho_times_C_rho ymm0
	vmovddup xmm0, xmm0 ; ymm0 = -|-|rho|rho
	vinsertf128 mask_rho_times_C_rho, mask_rho_times_C_rho, xmm0, 1 ; ymm0 = rho|rho|rho|rho
	
	%define mask_delta_t ymm1
	vmovddup delta_t, delta_t ; ymm1 = -|-|delta_t|delta_t
	vinsertf128 mask_delta_t, mask_delta_t, delta_t, 1 ; ymm1 = delta_t|delta_t|delta_t|delta_t
	
	%define i r8
	mov i, max_ij
	sub i, 4 ; i = max_ij - 4
	
	; primero calculamos los ultimos 4 elementos de la matriz
	shl i, 3
	vmovupd ymm2, [Tn + i]
	vmulpd ymm2, mask_rho_times_C_rho
	vdivpd ymm2, mask_delta_t
	vsubpd ymm2, [TIndAct + i] ; ymm2 = Tn * (- rho * C_rho) / delta_t - TIndAct
	vmovupd [B + i], ymm2
	shr i, 3
	
	; luego vemos donde comenzar el loop descendiente: debe ser el max k
	; tal que 4k <= max_ij - 1 - 4 (el ultimo 4k que no esta en los ultimos 4)
	dec i
	shr i, 2 ; i = k
	shl i, 5 ; i = 4k (*8 bytes)

	add i, 32 ; para contrarrestar lo que voy a restar despues
	loop_i:
		sub i, 32 ; saco los 4 elementos ya analizados (*8 bytes)
		
		vmovupd ymm2, [Tn + i]
		vmulpd ymm2, mask_rho_times_C_rho
		vdivpd ymm2, mask_delta_t
		vsubpd ymm2, [TIndAct + i] ; ymm2 = Tn * (- rho * C_rho) / delta_t - TIndAct
		vmovupd [B + i], ymm2
		
		jnz loop_i 
	
	pop rbp
ret
