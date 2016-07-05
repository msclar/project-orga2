section .data
    zeros: dq 0.0, 0.0

section .text
	global distanceBetweenVectors

distanceBetweenVectors:
	push rbp
	mov rbp, rsp
	
	%define actPhi  rdi
	%define nextPhi rsi
	%define length  rdx
	
	%define accum ymm0
	vxorpd accum, accum
	
	%define i r9
	mov i, length
	dec i
	shr i, 2
	shl i, 2 ; i = el mayor indice multiplo de 4 
	
	%define first_elem_not_analyzed r8
	mov first_elem_not_analyzed, i
	
	sub i, 4
	
	shl i, 3 ; pues cada double tiene 8 bytes
	shl first_elem_not_analyzed, 3
	shl length, 3

	add i, 32 ; para contrarrestar lo que voy a restar despues
	loop_i:
		sub i, 32 ; saco los 4 elementos ya analizados (*8 bytes)
		
		; accum += (actPhi[i] - nextPhi[i]) * (actPhi[i] - nextPhi[i]);
		vmovupd ymm1, [actPhi + i]
		vsubpd ymm1, [nextPhi + i] 
		vmulpd ymm1, ymm1
		vaddpd accum, ymm1

		jnz loop_i 
	
	
	%define unpacked_accum xmm1
	
	; desempaqueto el acumulador
	vextractf128 unpacked_accum, accum, 0x01 ; extraigo los 128 bits superiores, los inferiores estan en xmm0 (alli esta accum)
	addpd unpacked_accum, xmm0
	haddpd unpacked_accum, unpacked_accum
	
	; analizo los elementos [first_elem_not_analyzed, length)
	loop_end:
		movq xmm0, [actPhi + first_elem_not_analyzed]
		subsd xmm0, [nextPhi + first_elem_not_analyzed]
		mulsd xmm0, xmm0
		addsd unpacked_accum, xmm0
		
		add first_elem_not_analyzed, 8
		cmp first_elem_not_analyzed, length
		jl loop_end
	
	sqrtsd unpacked_accum, unpacked_accum
	movq xmm0, unpacked_accum
	
	pop rbp
ret
