section .data
    format: db "num: %d" , 10, 0

section .text
    global fillWithConstant

fillWithConstant:
    push rbp
    mov rbp, rsp
    
    %define matrix   rdi
    %define max_ij   rsi
    
    %define constant xmm0
    movddup constant, constant
    
    %define mask_constant ymm0
    vinsertf128 mask_constant, mask_constant, constant, 1

    %define i r8
    mov i, max_ij
    sub i, 4 ; i = max_ij - 4
    
    ; proceso los ultimos 4 elementos
    vmovupd [matrix + 8*i], mask_constant
    
    ; redondeo para abajo a un multiplo de 4 (idem updateB.asm)
    dec i
    shr i, 2
    shl i, 5

    add i, 32 ; para contrarrestar lo que voy a restar despues
    loop_i:
        sub i, 32 ; saco los 4 elementos ya analizados (*8 bytes)
        vmovupd [matrix + i], mask_constant
        jnz loop_i 
    
    pop rbp
ret
