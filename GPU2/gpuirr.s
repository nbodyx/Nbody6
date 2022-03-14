	.file	"gpuirr.avx.cpp"
	.text
	.p2align 4,,15
	.type	_ZL11gpuirr_openii._omp_fn.0, @function
_ZL11gpuirr_openii._omp_fn.0:
.LFB2353:
	.cfi_startproc
	ret
	.cfi_endproc
.LFE2353:
	.size	_ZL11gpuirr_openii._omp_fn.0, .-_ZL11gpuirr_openii._omp_fn.0
	.p2align 4,,15
	.type	_ZL11gpuirr_firriPdS_, @function
_ZL11gpuirr_firriPdS_:
.LFB1863:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movslq	%edi, %rdi
	leaq	(%rdi,%rdi,2), %rax
	salq	$5, %rax
	imulq	$3216, %rdi, %r9
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	movq	_ZL4ptcl(%rip), %r8
	vmovapd	_ZL8vec_tnow(%rip), %xmm15
	addq	%r8, %rax
	vmovddup	80(%rax), %xmm1
	vsubpd	%xmm1, %xmm15, %xmm1
	vbroadcastf128	48(%rax), %ymm2
	vbroadcastf128	64(%rax), %ymm3
	vcvtpd2psx	%xmm1, %xmm1
	vshufps	$0, %xmm1, %xmm1, %xmm0
	vshufps	$85, %xmm1, %xmm1, %xmm1
	vinsertf128	$0x1, %xmm1, %ymm0, %ymm0
	vmovaps	%ymm0, %ymm23
	vfmadd132ps	%ymm3, %ymm2, %ymm23
	vbroadcastf128	32(%rax), %ymm1
	vbroadcastf128	16(%rax), %ymm4
	vmovaps	.LC0(%rip), %ymm16
	vaddps	%ymm0, %ymm0, %ymm24
	vfmadd132ps	%ymm0, %ymm1, %ymm23
	addq	_ZL4list(%rip), %r9
	vbroadcastf32x4	(%rax), %ymm27
	movl	12(%r9), %ecx
	vshufps	$0, %ymm27, %ymm27, %ymm28
	vfmadd132ps	%ymm0, %ymm4, %ymm23
	vmulps	%ymm16, %ymm0, %ymm0
	vshufps	$85, %ymm27, %ymm27, %ymm29
	vmovaps	%ymm28, -32(%rsp)
	vmovaps	%ymm29, -64(%rsp)
	vshufps	$0, %ymm23, %ymm23, %ymm28
	vfmadd132ps	%ymm3, %ymm2, %ymm0
	vshufps	$85, %ymm23, %ymm23, %ymm29
	vshufps	$170, %ymm27, %ymm27, %ymm27
	vshufps	$170, %ymm23, %ymm23, %ymm23
	vfmadd132ps	%ymm0, %ymm1, %ymm24
	vshufps	$0, %ymm24, %ymm24, %ymm30
	vshufps	$85, %ymm24, %ymm24, %ymm31
	vshufps	$170, %ymm24, %ymm24, %ymm24
	testl	%ecx, %ecx
	jle	.L7
	leal	-1(%rcx), %eax
	leaq	16(%r9), %rdi
	movq	%rax, %r10
	leaq	20(%r9,%rax,4), %r11
	movq	%rdi, %rcx
	.p2align 4,,10
	.p2align 3
.L5:
	movslq	(%rcx), %rax
	addq	$4, %rcx
	leaq	(%rax,%rax,2), %rax
	salq	$5, %rax
	addq	%r8, %rax
	prefetcht0	(%rax)
	prefetcht0	64(%rax)
	cmpq	%rcx, %r11
	jne	.L5
	movl	%r10d, %ecx
	shrl	$3, %ecx
	vxorps	%xmm17, %xmm17, %xmm17
	salq	$5, %rcx
	vmovaps	.LC1(%rip), %ymm25
	leaq	48(%r9,%rcx), %r9
	vmovaps	%ymm17, %ymm18
	vmovaps	%ymm17, %ymm19
	vmovaps	%ymm17, %ymm20
	vmovaps	%ymm17, %ymm21
	vmovaps	%ymm17, %ymm22
	.p2align 4,,10
	.p2align 3
.L6:
	movslq	(%rdi), %rcx
	movslq	4(%rdi), %rax
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rcx
	addq	%r8, %rcx
	leaq	(%rax,%rax,2), %rax
	vmovapd	80(%rcx), %xmm26
	salq	$5, %rax
	addq	%r8, %rax
	vmovhpd	80(%rax), %xmm26, %xmm9
	vmovaps	16(%rcx), %xmm6
	vmovaps	32(%rcx), %xmm4
	vsubpd	%xmm9, %xmm15, %xmm9
	vinsertf128	$0x1, 16(%rax), %ymm6, %ymm2
	vmovaps	48(%rcx), %xmm3
	vinsertf128	$0x1, 32(%rax), %ymm4, %ymm6
	vmovaps	64(%rcx), %xmm4
	vinsertf128	$0x1, 48(%rax), %ymm3, %ymm12
	vcvtpd2psx	%xmm9, %xmm9
	vshufps	$0, %xmm9, %xmm9, %xmm0
	vshufps	$85, %xmm9, %xmm9, %xmm9
	vinsertf128	$0x1, %xmm9, %ymm0, %ymm0
	vinsertf128	$0x1, 64(%rax), %ymm4, %ymm9
	vmovaps	%ymm0, %ymm8
	vfmadd132ps	%ymm9, %ymm12, %ymm8
	vmovaps	(%rcx), %xmm11
	movslq	8(%rdi), %rcx
	vinsertf128	$0x1, (%rax), %ymm11, %ymm11
	leaq	(%rcx,%rcx,2), %rcx
	vfmadd132ps	%ymm0, %ymm6, %ymm8
	movslq	12(%rdi), %rax
	salq	$5, %rcx
	addq	%r8, %rcx
	leaq	(%rax,%rax,2), %rax
	vfmadd132ps	%ymm0, %ymm2, %ymm8
	vaddps	%ymm0, %ymm0, %ymm2
	vmulps	%ymm16, %ymm0, %ymm0
	vmovapd	80(%rcx), %xmm26
	salq	$5, %rax
	addq	%r8, %rax
	vmovaps	64(%rcx), %xmm3
	vfmadd132ps	%ymm9, %ymm12, %ymm0
	vmovaps	48(%rcx), %xmm7
	vmovaps	(%rcx), %xmm10
	vinsertf128	$0x1, 48(%rax), %ymm7, %ymm7
	vinsertf128	$0x1, (%rax), %ymm10, %ymm10
	vfmadd132ps	%ymm0, %ymm6, %ymm2
	vmovhpd	80(%rax), %xmm26, %xmm6
	vsubpd	%xmm6, %xmm15, %xmm6
	addq	$32, %rdi
	vmovaps	%ymm2, %ymm9
	vcvtpd2psx	%xmm6, %xmm6
	vmovaps	16(%rcx), %xmm2
	vshufps	$0, %xmm6, %xmm6, %xmm0
	vshufps	$85, %xmm6, %xmm6, %xmm6
	vinsertf128	$0x1, %xmm6, %ymm0, %ymm0
	vmovaps	32(%rcx), %xmm6
	vinsertf128	$0x1, 16(%rax), %ymm2, %ymm1
	vinsertf128	$0x1, 32(%rax), %ymm6, %ymm2
	vinsertf128	$0x1, 64(%rax), %ymm3, %ymm6
	vmovaps	%ymm0, %ymm4
	vfmadd132ps	%ymm6, %ymm7, %ymm4
	movslq	-16(%rdi), %rcx
	movslq	-12(%rdi), %rax
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rcx
	vfmadd132ps	%ymm0, %ymm2, %ymm4
	addq	%r8, %rcx
	leaq	(%rax,%rax,2), %rax
	vmovapd	80(%rcx), %xmm26
	salq	$5, %rax
	vfmadd132ps	%ymm0, %ymm1, %ymm4
	vaddps	%ymm0, %ymm0, %ymm1
	vmulps	%ymm16, %ymm0, %ymm0
	addq	%r8, %rax
	vfmadd132ps	%ymm6, %ymm7, %ymm0
	vmovaps	%ymm1, %ymm6
	vmovhpd	80(%rax), %xmm26, %xmm1
	vsubpd	%xmm1, %xmm15, %xmm1
	vfmadd132ps	%ymm0, %ymm2, %ymm6
	vmovaps	16(%rcx), %xmm2
	vmovaps	32(%rcx), %xmm7
	vmovaps	48(%rcx), %xmm3
	vinsertf128	$0x1, 32(%rax), %ymm7, %ymm5
	vinsertf128	$0x1, 48(%rax), %ymm3, %ymm7
	vmovaps	64(%rcx), %xmm3
	vcvtpd2psx	%xmm1, %xmm1
	vinsertf128	$0x1, 64(%rax), %ymm3, %ymm12
	vshufps	$0, %xmm1, %xmm1, %xmm0
	vshufps	$85, %xmm1, %xmm1, %xmm1
	vinsertf128	$0x1, %xmm1, %ymm0, %ymm0
	vmovaps	%ymm0, %ymm1
	vfmadd132ps	%ymm12, %ymm7, %ymm1
	vinsertf128	$0x1, 16(%rax), %ymm2, %ymm2
	vmovaps	(%rcx), %xmm3
	movslq	-8(%rdi), %rcx
	vinsertf128	$0x1, (%rax), %ymm3, %ymm3
	vfmadd132ps	%ymm0, %ymm5, %ymm1
	movslq	-4(%rdi), %rax
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rcx
	addq	%r8, %rcx
	vfmadd132ps	%ymm0, %ymm2, %ymm1
	vaddps	%ymm0, %ymm0, %ymm2
	vmulps	%ymm16, %ymm0, %ymm0
	leaq	(%rax,%rax,2), %rax
	vmovapd	80(%rcx), %xmm26
	salq	$5, %rax
	addq	%r8, %rax
	vfmadd132ps	%ymm12, %ymm7, %ymm0
	vmovaps	32(%rcx), %xmm7
	vinsertf128	$0x1, 32(%rax), %ymm7, %ymm12
	vmovaps	48(%rcx), %xmm7
	vfmadd132ps	%ymm0, %ymm5, %ymm2
	vmovhpd	80(%rax), %xmm26, %xmm5
	vsubpd	%xmm5, %xmm15, %xmm5
	vinsertf128	$0x1, 48(%rax), %ymm7, %ymm13
	vmovaps	64(%rcx), %xmm7
	vcvtpd2psx	%xmm5, %xmm5
	vinsertf128	$0x1, 64(%rax), %ymm7, %ymm14
	vshufps	$0, %xmm5, %xmm5, %xmm0
	vshufps	$85, %xmm5, %xmm5, %xmm5
	vinsertf128	$0x1, %xmm5, %ymm0, %ymm0
	vmovaps	%ymm0, %ymm7
	vfmadd132ps	%ymm14, %ymm13, %ymm7
	vmovaps	16(%rcx), %xmm26
	vmovaps	(%rcx), %xmm5
	vinsertf32x4	$0x1, 16(%rax), %ymm26, %ymm26
	vinsertf128	$0x1, (%rax), %ymm5, %ymm5
	vfmadd132ps	%ymm0, %ymm12, %ymm7
	vfmadd132ps	%ymm0, %ymm26, %ymm7
	vaddps	%ymm0, %ymm0, %ymm26
	vmulps	%ymm16, %ymm0, %ymm0
	vfmadd132ps	%ymm14, %ymm13, %ymm0
	vunpcklps	%ymm3, %ymm11, %ymm14
	vunpcklps	%ymm5, %ymm10, %ymm13
	vunpckhps	%ymm3, %ymm11, %ymm3
	vunpckhps	%ymm5, %ymm10, %ymm10
	vfmadd231ps	%ymm0, %ymm26, %ymm12
	vunpcklps	%ymm1, %ymm8, %ymm5
	vunpcklps	%ymm13, %ymm14, %ymm11
	vunpckhps	%ymm13, %ymm14, %ymm13
	vunpcklps	%ymm10, %ymm3, %ymm14
	vunpckhps	%ymm10, %ymm3, %ymm3
	vunpcklps	%ymm7, %ymm4, %ymm10
	vunpcklps	%ymm10, %ymm5, %ymm0
	vunpckhps	%ymm10, %ymm5, %ymm10
	vsubps	%ymm29, %ymm10, %ymm10
	vsubps	-64(%rsp), %ymm13, %ymm13
	vunpckhps	%ymm7, %ymm4, %ymm4
	vunpckhps	%ymm1, %ymm8, %ymm1
	vunpcklps	%ymm12, %ymm6, %ymm5
	vunpcklps	%ymm2, %ymm9, %ymm7
	vunpcklps	%ymm4, %ymm1, %ymm1
	vunpcklps	%ymm5, %ymm7, %ymm8
	vsubps	%ymm23, %ymm1, %ymm1
	vunpckhps	%ymm5, %ymm7, %ymm7
	vsubps	%ymm27, %ymm14, %ymm4
	vaddps	%ymm10, %ymm13, %ymm5
	vsubps	%ymm28, %ymm0, %ymm0
	vsubps	-32(%rsp), %ymm11, %ymm11
	vaddps	%ymm1, %ymm4, %ymm4
	vmulps	%ymm5, %ymm5, %ymm1
	vaddps	%ymm0, %ymm11, %ymm11
	vsubps	%ymm31, %ymm7, %ymm7
	vsubps	%ymm30, %ymm8, %ymm8
	vunpckhps	%ymm2, %ymm9, %ymm2
	vfmadd231ps	%ymm11, %ymm11, %ymm1
	vunpckhps	%ymm12, %ymm6, %ymm12
	vunpcklps	%ymm12, %ymm2, %ymm12
	vsubps	%ymm24, %ymm12, %ymm12
	vfmadd231ps	%ymm4, %ymm4, %ymm1
	vrsqrt14ps	%ymm1, %ymm0
	vmulps	%ymm0, %ymm1, %ymm1
	vfmadd132ps	%ymm0, %ymm25, %ymm1
	vmulps	.LC2(%rip), %ymm0, %ymm0
	vmulps	%ymm0, %ymm1, %ymm0
	vmulps	%ymm7, %ymm5, %ymm1
	vmulps	%ymm0, %ymm0, %ymm10
	vfmadd231ps	%ymm8, %ymm11, %ymm1
	vmulps	%ymm3, %ymm0, %ymm3
	vmulps	%ymm25, %ymm10, %ymm2
	vfmadd231ps	%ymm12, %ymm4, %ymm1
	vmulps	%ymm10, %ymm3, %ymm3
	vmulps	%ymm2, %ymm1, %ymm1
	vfmadd231ps	%ymm3, %ymm11, %ymm22
	vfmadd231ps	%ymm3, %ymm5, %ymm21
	vfmadd231ps	%ymm3, %ymm4, %ymm20
	vfmadd132ps	%ymm1, %ymm8, %ymm11
	vfmadd132ps	%ymm1, %ymm7, %ymm5
	vfmadd132ps	%ymm1, %ymm12, %ymm4
	vfmadd231ps	%ymm11, %ymm3, %ymm19
	vfmadd231ps	%ymm5, %ymm3, %ymm18
	vfmadd231ps	%ymm4, %ymm3, %ymm17
	cmpq	%rdi, %r9
	jne	.L6
.L4:
	vcvtps2pd	%xmm22, %ymm1
	vextractf32x4	$0x1, %ymm22, %xmm22
	vcvtps2pd	%xmm22, %ymm22
	vaddpd	%ymm22, %ymm1, %ymm1
	vhaddpd	%ymm1, %ymm1, %ymm1
	vmovapd	%xmm1, %xmm0
	vextractf64x2	$0x1, %ymm1, %xmm1
	vaddpd	%xmm1, %xmm0, %xmm1
	vcvtps2pd	%xmm21, %ymm0
	vextractf32x4	$0x1, %ymm21, %xmm21
	vcvtps2pd	%xmm21, %ymm21
	vaddpd	%ymm21, %ymm0, %ymm0
	vmovlpd	%xmm1, (%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm21
	vextractf64x2	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm21, %xmm21
	vcvtps2pd	%xmm20, %ymm0
	vextractf32x4	$0x1, %ymm20, %xmm20
	vcvtps2pd	%xmm20, %ymm20
	vaddpd	%ymm20, %ymm0, %ymm0
	vmovlpd	%xmm21, 8(%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm20
	vextractf64x2	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm20, %xmm20
	vcvtps2pd	%xmm19, %ymm0
	vextractf32x4	$0x1, %ymm19, %xmm19
	vcvtps2pd	%xmm19, %ymm19
	vaddpd	%ymm19, %ymm0, %ymm0
	vmovlpd	%xmm20, 16(%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm19
	vextractf64x2	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm19, %xmm19
	vcvtps2pd	%xmm18, %ymm0
	vextractf32x4	$0x1, %ymm18, %xmm18
	vcvtps2pd	%xmm18, %ymm18
	vaddpd	%ymm18, %ymm0, %ymm0
	vmovlpd	%xmm19, (%rdx)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm18
	vextractf64x2	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm18, %xmm18
	vcvtps2pd	%xmm17, %ymm0
	vextractf32x4	$0x1, %ymm17, %xmm17
	vcvtps2pd	%xmm17, %ymm17
	vaddpd	%ymm17, %ymm0, %ymm0
	vmovlpd	%xmm18, 8(%rdx)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm2
	vextractf64x2	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm2, %xmm0
	vmovlpd	%xmm0, 16(%rdx)
	vzeroupper
	leave
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L7:
	.cfi_restore_state
	vxorps	%xmm17, %xmm17, %xmm17
	vmovaps	%ymm17, %ymm18
	vmovaps	%ymm17, %ymm19
	vmovaps	%ymm17, %ymm20
	vmovaps	%ymm17, %ymm21
	vmovaps	%ymm17, %ymm22
	jmp	.L4
	.cfi_endproc
.LFE1863:
	.size	_ZL11gpuirr_firriPdS_, .-_ZL11gpuirr_firriPdS_
	.p2align 4,,15
	.type	_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.1, @function
_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.1:
.LFB2354:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movl	$1, %edx
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	xorl	%ebx, %ebx
	subq	$56, %rsp
	.cfi_def_cfa_offset 112
	movq	16(%rdi), %rcx
	movslq	24(%rdi), %rsi
	movq	%rcx, 8(%rsp)
	movq	8(%rdi), %rcx
	movq	%rdi, 24(%rsp)
	movq	%rcx, 16(%rsp)
	movq	(%rdi), %r14
	leaq	40(%rsp), %r9
	leaq	32(%rsp), %r8
	movl	$1, %ecx
	xorl	%edi, %edi
	call	GOMP_loop_guided_start
	testb	%al, %al
	je	.L12
	movq	_ZL4list(%rip), %r12
	.p2align 4,,10
	.p2align 3
.L14:
	movl	40(%rsp), %eax
	movslq	32(%rsp), %rbp
	movl	%eax, 4(%rsp)
	leaq	0(%rbp,%rbp,2), %rcx
	movq	16(%rsp), %rax
	salq	$3, %rcx
	leaq	(%rax,%rcx), %r13
	movq	8(%rsp), %rax
	leaq	(%rax,%rcx), %r15
	.p2align 4,,10
	.p2align 3
.L13:
	movl	(%r14,%rbp,4), %eax
	movq	%r15, %rdx
	movq	%r13, %rsi
	leal	-1(%rax), %edi
	call	_ZL11gpuirr_firriPdS_
	movslq	(%r14,%rbp,4), %rdx
	incq	%rbp
	imulq	$3216, %rdx, %rdx
	addq	$24, %r13
	addq	$24, %r15
	addl	-3204(%r12,%rdx), %ebx
	cmpl	%ebp, 4(%rsp)
	jg	.L13
	leaq	40(%rsp), %rsi
	leaq	32(%rsp), %rdi
	call	GOMP_loop_guided_next
	testb	%al, %al
	jne	.L14
.L12:
	call	GOMP_loop_end_nowait
	movq	24(%rsp), %rax
	lock addl	%ebx, 28(%rax)
	addq	$56, %rsp
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2354:
	.size	_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.1, .-_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.1
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"gpuirr: it is already open\n"
.LC4:
	.string	"irrlib/gpuirr.avx.cpp"
.LC5:
	.string	"lmax <= 1 + NBlist::NB_MAX"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC6:
	.string	"**************************** \n"
	.align 8
.LC7:
	.string	"Opening GPUIRR lib. AVX ver. \n"
	.section	.rodata.str1.1
.LC8:
	.string	" nmax = %d, lmax = %d\n"
	.section	.rodata.str1.8
	.align 8
.LC11:
	.string	"0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist))"
	.align 8
.LC12:
	.string	"0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle))"
	.text
	.p2align 4,,15
	.globl	gpuirr_open_
	.type	gpuirr_open_, @function
gpuirr_open_:
.LFB1865:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$24, %rsp
	.cfi_def_cfa_offset 48
	cmpb	$0, _ZL7is_open(%rip)
	jne	.L28
	movl	(%rsi), %ebx
	cmpl	$801, %ebx
	jg	.L29
	movq	stderr(%rip), %rcx
	movl	(%rdi), %ebp
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC7, %edi
	call	fwrite
	movq	stderr(%rip), %rdi
	movl	%ebx, %ecx
	movl	%ebp, %edx
	movl	$.LC8, %esi
	xorl	%eax, %eax
	call	fprintf
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	leal	1(%rbp), %eax
	cltq
	leaq	(%rax,%rax,2), %rbx
	salq	$5, %rbx
	movq	%rbx, %rdx
	movl	$64, %esi
	leaq	8(%rsp), %rdi
	call	posix_memalign
	testl	%eax, %eax
	jne	.L24
	movq	8(%rsp), %rcx
	movq	%rbx, %rdx
	movl	$255, %esi
	movq	%rcx, %rdi
	call	memset
	movq	%rax, _ZL4ptcl(%rip)
	vmovaps	.LC9(%rip), %xmm0
	leaq	-96(%rax,%rbx), %rax
	movslq	%ebp, %rbx
	imulq	$3216, %rbx, %rbx
	vmovaps	%xmm0, (%rax)
	vxorps	%xmm0, %xmm0, %xmm0
	vmovaps	%xmm0, 16(%rax)
	vmovaps	%xmm0, 32(%rax)
	vmovaps	%xmm0, 48(%rax)
	vmovaps	%xmm0, 64(%rax)
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovaps	%xmm0, 80(%rax)
	movq	%rbx, %rdx
	movl	$64, %esi
	movq	%rsp, %rdi
	call	posix_memalign
	testl	%eax, %eax
	jne	.L25
	movq	(%rsp), %rcx
	movq	%rbx, %rdx
	movq	%rcx, %rdi
	movl	$255, %esi
	call	memset
	xorl	%ecx, %ecx
	xorl	%edx, %edx
	xorl	%esi, %esi
	movl	$_ZL11gpuirr_openii._omp_fn.0, %edi
	movl	%ebp, _ZL4nmax(%rip)
	movq	%rax, _ZL4list(%rip)
	call	GOMP_parallel
	movq	$0x000000000, _ZL9time_grav(%rip)
	movq	$0, _ZL9num_steps(%rip)
	movq	$0, _ZL9num_fcall(%rip)
	movq	$0, _ZL9num_inter(%rip)
	movb	$1, _ZL7is_open(%rip)
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L28:
	.cfi_restore_state
	movq	stderr(%rip), %rcx
	movl	$27, %edx
	movl	$1, %esi
	movl	$.LC3, %edi
	call	fwrite
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L24:
	.cfi_restore_state
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$370, %edx
	movl	$.LC4, %esi
	movl	$.LC12, %edi
	call	__assert_fail
.L25:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$378, %edx
	movl	$.LC4, %esi
	movl	$.LC11, %edi
	call	__assert_fail
.L29:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$360, %edx
	movl	$.LC4, %esi
	movl	$.LC5, %edi
	call	__assert_fail
	.cfi_endproc
.LFE1865:
	.size	gpuirr_open_, .-gpuirr_open_
	.section	.rodata.str1.1
.LC13:
	.string	"gpuirr: it is already close\n"
	.section	.rodata.str1.8
	.align 8
.LC15:
	.string	"Closing GPUIRR lib. AVX ver. \n"
	.section	.rodata.str1.1
.LC16:
	.string	"time grav  : %f sec\n"
.LC19:
	.string	"perf grav  : %f Gflops\n"
.LC20:
	.string	"perf grav  : %f usec\n"
.LC21:
	.string	"<#NB>      : %f \n"
	.text
	.p2align 4,,15
	.globl	gpuirr_close_
	.type	gpuirr_close_, @function
gpuirr_close_:
.LFB1866:
	.cfi_startproc
	cmpb	$0, _ZL7is_open(%rip)
	je	.L40
	subq	$40, %rsp
	.cfi_def_cfa_offset 48
	movq	_ZL4ptcl(%rip), %rdi
	call	free
	movq	_ZL4list(%rip), %rdi
	movq	$0, _ZL4ptcl(%rip)
	call	free
	movq	_ZL9num_inter(%rip), %rax
	movq	$0, _ZL4list(%rip)
	testq	%rax, %rax
	js	.L32
	vxorpd	%xmm1, %xmm1, %xmm1
	vcvtsi2sdq	%rax, %xmm1, %xmm1
.L33:
	movq	_ZL9num_fcall(%rip), %rax
	vmovsd	_ZL9time_grav(%rip), %xmm4
	testq	%rax, %rax
	js	.L34
	vxorpd	%xmm2, %xmm2, %xmm2
	vcvtsi2sdq	%rax, %xmm2, %xmm2
.L35:
	vdivsd	%xmm2, %xmm4, %xmm2
	movq	_ZL9num_steps(%rip), %rax
	vmulsd	.LC14(%rip), %xmm2, %xmm2
	testq	%rax, %rax
	js	.L36
	vxorpd	%xmm3, %xmm3, %xmm3
	vcvtsi2sdq	%rax, %xmm3, %xmm3
.L37:
	vdivsd	%xmm3, %xmm1, %xmm3
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	vmovsd	%xmm2, 24(%rsp)
	vmovsd	%xmm4, 16(%rsp)
	vmovsd	%xmm1, 8(%rsp)
	vmovsd	%xmm3, (%rsp)
	call	fwrite
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC15, %edi
	call	fwrite
	vmovsd	_ZL9time_grav(%rip), %xmm0
	movq	stderr(%rip), %rdi
	movl	$.LC16, %esi
	movl	$1, %eax
	call	fprintf
	movq	stderr(%rip), %rsi
	movl	$10, %edi
	call	fputc
	vmovsd	8(%rsp), %xmm1
	vmovsd	16(%rsp), %xmm4
	vmulsd	.LC17(%rip), %xmm1, %xmm0
	movq	stderr(%rip), %rdi
	movl	$.LC19, %esi
	movl	$1, %eax
	vmulsd	.LC18(%rip), %xmm0, %xmm0
	vdivsd	%xmm4, %xmm0, %xmm0
	call	fprintf
	vmovsd	24(%rsp), %xmm2
	movq	stderr(%rip), %rdi
	vmovapd	%xmm2, %xmm0
	movl	$.LC20, %esi
	movl	$1, %eax
	call	fprintf
	vmovsd	(%rsp), %xmm3
	movq	stderr(%rip), %rdi
	vmovapd	%xmm3, %xmm0
	movl	$.LC21, %esi
	movl	$1, %eax
	call	fprintf
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	movb	$0, _ZL7is_open(%rip)
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L32:
	.cfi_restore_state
	movq	%rax, %rdx
	shrq	%rdx
	andl	$1, %eax
	orq	%rax, %rdx
	vxorpd	%xmm1, %xmm1, %xmm1
	vcvtsi2sdq	%rdx, %xmm1, %xmm1
	vaddsd	%xmm1, %xmm1, %xmm1
	jmp	.L33
	.p2align 4,,10
	.p2align 3
.L40:
	.cfi_def_cfa_offset 8
	movq	stderr(%rip), %rcx
	movl	$28, %edx
	movl	$1, %esi
	movl	$.LC13, %edi
	jmp	fwrite
	.p2align 4,,10
	.p2align 3
.L36:
	.cfi_def_cfa_offset 48
	movq	%rax, %rdx
	shrq	%rdx
	andl	$1, %eax
	orq	%rax, %rdx
	vxorpd	%xmm3, %xmm3, %xmm3
	vcvtsi2sdq	%rdx, %xmm3, %xmm3
	vaddsd	%xmm3, %xmm3, %xmm3
	jmp	.L37
	.p2align 4,,10
	.p2align 3
.L34:
	movq	%rax, %rdx
	shrq	%rdx
	andl	$1, %eax
	orq	%rax, %rdx
	vxorpd	%xmm2, %xmm2, %xmm2
	vcvtsi2sdq	%rdx, %xmm2, %xmm2
	vaddsd	%xmm2, %xmm2, %xmm2
	jmp	.L35
	.cfi_endproc
.LFE1866:
	.size	gpuirr_close_, .-gpuirr_close_
	.p2align 4,,15
	.globl	gpuirr_set_jp_
	.type	gpuirr_set_jp_, @function
gpuirr_set_jp_:
.LFB1867:
	.cfi_startproc
	vmovsd	16(%rsi), %xmm1
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	vmovupd	(%rsi), %xmm0
	vmovhpd	(%r9), %xmm1, %xmm1
	vinsertf128	$0x1, %xmm1, %ymm0, %ymm0
	vmovupd	(%rdx), %xmm2
	vcvtpd2psy	%ymm0, %xmm4
	vcvtps2pd	%xmm4, %ymm1
	vsubpd	%ymm1, %ymm0, %ymm0
	vmovq	16(%rdx), %xmm1
	movl	(%rdi), %eax
	vmovupd	(%rcx), %xmm3
	vinsertf128	$0x1, %xmm1, %ymm2, %ymm2
	vmovq	16(%rcx), %xmm1
	decl	%eax
	vmovq	16(%r8), %xmm5
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	vinsertf128	$0x1, %xmm1, %ymm3, %ymm3
	cltq
	vmovupd	(%r8), %xmm1
	movq	16(%rbp), %rdx
	leaq	(%rax,%rax,2), %rax
	vinsertf128	$0x1, %xmm5, %ymm1, %ymm1
	salq	$5, %rax
	addq	_ZL4ptcl(%rip), %rax
	vcvtpd2psy	%ymm3, %xmm3
	vcvtpd2psy	%ymm0, %xmm0
	vcvtpd2psy	%ymm2, %xmm2
	vinsertf128	$0x1, %xmm0, %ymm4, %ymm0
	vinsertf128	$0x1, %xmm3, %ymm2, %ymm2
	vcvtpd2psy	%ymm1, %xmm1
	vmovddup	(%rdx), %xmm3
	vinsertf128	$0x1, %xmm3, %ymm1, %ymm1
	vmovaps	%ymm0, (%rax)
	vmovaps	%ymm2, 32(%rax)
	vmovaps	%ymm1, 64(%rax)
	vzeroupper
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE1867:
	.size	gpuirr_set_jp_, .-gpuirr_set_jp_
	.section	.rodata.str1.1
.LC22:
	.string	"nnb <= NBlist::NB_MAX"
	.text
	.p2align 4,,15
	.globl	gpuirr_set_list_
	.type	gpuirr_set_list_, @function
gpuirr_set_list_:
.LFB1868:
	.cfi_startproc
	movl	(%rdi), %edx
	movl	(%rsi), %eax
	decl	%edx
	cmpl	$800, %eax
	jg	.L60
	movslq	%edx, %rdx
	imulq	$3216, %rdx, %rdx
	leal	-1(%rax), %r10d
	addq	_ZL4list(%rip), %rdx
	movl	%eax, 12(%rdx)
	leaq	16(%rdx), %rdi
	testl	%eax, %eax
	jle	.L45
	leal	-1(%rax), %r10d
	movl	%r10d, %r8d
	shrl	$3, %r8d
	salq	$5, %r8
	vmovdqa64	.LC23(%rip), %xmm2
	leaq	4(%rsi), %rcx
	leaq	32(%rdx), %r9
	leaq	36(%rsi,%r8), %rsi
	.p2align 4,,10
	.p2align 3
.L46:
	vmovdqu8	16(%rcx), %xmm0
	vmovdqu32	(%rcx), %xmm3
	vpsubd	%xmm2, %xmm0, %xmm0
	vpsubd	%xmm2, %xmm3, %xmm1
	addq	$32, %rcx
	vmovntdq	%xmm1, -16(%r9)
	vmovntdq	%xmm0, (%r9)
	addq	$32, %r9
	cmpq	%rsi, %rcx
	jne	.L46
.L45:
	testl	%r10d, %r10d
	leal	6(%rax), %ecx
	cmovns	%r10d, %ecx
	movl	_ZL4nmax(%rip), %esi
	sarl	$3, %ecx
	leal	8(,%rcx,8), %ecx
	cmpl	%ecx, %eax
	jge	.L56
	movl	%eax, %r8d
	notl	%r8d
	movl	%ecx, %r9d
	addl	%ecx, %r8d
	subl	%eax, %r9d
	cmpl	$6, %r8d
	jbe	.L48
	movslq	%eax, %r8
	leaq	16(%rdx,%r8,4), %r8
	movl	%r9d, %edx
	shrl	$3, %edx
	salq	$5, %rdx
	vpbroadcastd	%esi, %ymm0
	addq	%r8, %rdx
	.p2align 4,,10
	.p2align 3
.L49:
	vmovdqu32	%ymm0, (%r8)
	addq	$32, %r8
	cmpq	%rdx, %r8
	jne	.L49
	movl	%r9d, %edx
	andl	$-8, %edx
	addl	%edx, %eax
	cmpl	%edx, %r9d
	je	.L55
	vzeroupper
.L48:
	movslq	%eax, %rdx
	movl	%esi, (%rdi,%rdx,4)
	leal	1(%rax), %edx
	cmpl	%edx, %ecx
	jle	.L56
	movslq	%edx, %rdx
	movl	%esi, (%rdi,%rdx,4)
	leal	2(%rax), %edx
	cmpl	%edx, %ecx
	jle	.L56
	movslq	%edx, %rdx
	movl	%esi, (%rdi,%rdx,4)
	leal	3(%rax), %edx
	cmpl	%edx, %ecx
	jle	.L56
	movslq	%edx, %rdx
	movl	%esi, (%rdi,%rdx,4)
	leal	4(%rax), %edx
	cmpl	%edx, %ecx
	jle	.L56
	movslq	%edx, %rdx
	movl	%esi, (%rdi,%rdx,4)
	leal	5(%rax), %edx
	cmpl	%edx, %ecx
	jle	.L56
	movslq	%edx, %rdx
	addl	$6, %eax
	movl	%esi, (%rdi,%rdx,4)
	cmpl	%eax, %ecx
	jle	.L56
	cltq
	movl	%esi, (%rdi,%rax,4)
	ret
	.p2align 4,,10
	.p2align 3
.L55:
	vzeroupper
.L56:
	ret
.L60:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movl	$_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, %ecx
	movl	$435, %edx
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	$.LC4, %esi
	andq	$-32, %rsp
	movl	$.LC22, %edi
	call	__assert_fail
	.cfi_endproc
.LFE1868:
	.size	gpuirr_set_list_, .-gpuirr_set_list_
	.p2align 4,,15
	.globl	gpuirr_pred_all_
	.type	gpuirr_pred_all_, @function
gpuirr_pred_all_:
.LFB1869:
	.cfi_startproc
	vmovddup	(%rdx), %xmm0
	vmovaps	%xmm0, _ZL8vec_tnow(%rip)
	ret
	.cfi_endproc
.LFE1869:
	.size	gpuirr_pred_all_, .-gpuirr_pred_all_
	.p2align 4,,15
	.globl	gpuirr_pred_act_
	.type	gpuirr_pred_act_, @function
gpuirr_pred_act_:
.LFB2359:
	.cfi_startproc
	vmovddup	(%rdx), %xmm0
	vmovaps	%xmm0, _ZL8vec_tnow(%rip)
	ret
	.cfi_endproc
.LFE2359:
	.size	gpuirr_pred_act_, .-gpuirr_pred_act_
	.p2align 4,,15
	.globl	gpuirr_firr_vec_
	.type	gpuirr_firr_vec_, @function
gpuirr_firr_vec_:
.LFB1871:
	.cfi_startproc
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	movq	%rcx, %r13
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	movq	%rdx, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	movq	%rsi, %rbp
	xorl	%esi, %esi
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$56, %rsp
	.cfi_def_cfa_offset 96
	movslq	(%rdi), %rbx
	leaq	16(%rsp), %rdi
	call	gettimeofday
	vxorpd	%xmm1, %xmm1, %xmm1
	vxorpd	%xmm0, %xmm0, %xmm0
	vcvtsi2sdq	16(%rsp), %xmm0, %xmm0
	vcvtsi2sdq	24(%rsp), %xmm1, %xmm1
	leaq	16(%rsp), %rsi
	xorl	%ecx, %ecx
	xorl	%edx, %edx
	vfmadd132sd	.LC24(%rip), %xmm0, %xmm1
	movl	$_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.1, %edi
	movq	%r13, 32(%rsp)
	movq	%r12, 24(%rsp)
	movq	%rbp, 16(%rsp)
	vmovsd	%xmm1, 8(%rsp)
	movl	%ebx, 40(%rsp)
	movl	$0, 44(%rsp)
	call	GOMP_parallel
	movslq	44(%rsp), %rax
	leaq	16(%rsp), %rdi
	xorl	%esi, %esi
	addq	%rax, _ZL9num_inter(%rip)
	call	gettimeofday
	vxorpd	%xmm0, %xmm0, %xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	vcvtsi2sdq	24(%rsp), %xmm0, %xmm0
	vcvtsi2sdq	16(%rsp), %xmm1, %xmm1
	addq	%rbx, _ZL9num_steps(%rip)
	incq	_ZL9num_fcall(%rip)
	vfmadd132sd	.LC24(%rip), %xmm1, %xmm0
	vsubsd	8(%rsp), %xmm0, %xmm0
	vaddsd	_ZL9time_grav(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, _ZL9time_grav(%rip)
	addq	$56, %rsp
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1871:
	.size	gpuirr_firr_vec_, .-gpuirr_firr_vec_
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.type	_GLOBAL__sub_I_gpuirr_open_, @function
_GLOBAL__sub_I_gpuirr_open_:
.LFB2352:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	jmp	__cxa_atexit
	.cfi_endproc
.LFE2352:
	.size	_GLOBAL__sub_I_gpuirr_open_, .-_GLOBAL__sub_I_gpuirr_open_
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I_gpuirr_open_
	.section	.rodata
	.align 32
	.type	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, 43
_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_set_list(int, int, const int*)"
	.align 16
	.type	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, 27
_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_open(int, int)"
	.local	_ZL8vec_tnow
	.comm	_ZL8vec_tnow,16,16
	.local	_ZL9num_steps
	.comm	_ZL9num_steps,8,8
	.local	_ZL9num_fcall
	.comm	_ZL9num_fcall,8,8
	.local	_ZL9num_inter
	.comm	_ZL9num_inter,8,8
	.local	_ZL9time_grav
	.comm	_ZL9time_grav,8,8
	.local	_ZL4nmax
	.comm	_ZL4nmax,4,4
	.local	_ZL4list
	.comm	_ZL4list,8,8
	.local	_ZL4ptcl
	.comm	_ZL4ptcl,8,8
	.local	_ZL7is_open
	.comm	_ZL7is_open,1,1
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC0:
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.align 32
.LC1:
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.align 32
.LC2:
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC9:
	.long	1132396544
	.long	1132396544
	.long	1132396544
	.long	0
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC14:
	.long	0
	.long	1093567616
	.align 8
.LC17:
	.long	0
	.long	1078853632
	.align 8
.LC18:
	.long	3894859413
	.long	1041313291
	.section	.rodata.cst16
	.align 16
.LC23:
	.long	1
	.long	1
	.long	1
	.long	1
	.section	.rodata.cst8
	.align 8
.LC24:
	.long	2696277389
	.long	1051772663
	.hidden	__dso_handle
	.ident	"GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-4)"
	.section	.note.GNU-stack,"",@progbits
