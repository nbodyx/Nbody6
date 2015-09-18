	.file	"gpuirr.avx.cpp"
	.text
	.p2align 4,,15
.globl gpuirr_pred_all_
	.type	gpuirr_pred_all_, @function
gpuirr_pred_all_:
.LFB2403:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	movq	(%rdx), %rax
	movq	%rax, _ZL8vec_tnow(%rip)
	movq	%rax, _ZL8vec_tnow+8(%rip)
	ret
	.cfi_endproc
.LFE2403:
	.size	gpuirr_pred_all_, .-gpuirr_pred_all_
	.p2align 4,,15
.globl gpuirr_pred_act_
	.type	gpuirr_pred_act_, @function
gpuirr_pred_act_:
.LFB2404:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	movq	(%rdx), %rax
	movq	%rax, _ZL8vec_tnow(%rip)
	movq	%rax, _ZL8vec_tnow+8(%rip)
	ret
	.cfi_endproc
.LFE2404:
	.size	gpuirr_pred_act_, .-gpuirr_pred_act_
	.p2align 4,,15
	.type	_GLOBAL__I_gpuirr_open_, @function
_GLOBAL__I_gpuirr_open_:
.LFB2510:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
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
.LFE2510:
	.size	_GLOBAL__I_gpuirr_open_, .-_GLOBAL__I_gpuirr_open_
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I_gpuirr_open_
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"irrlib/gpuirr.avx.cpp"
.LC1:
	.string	"nnb <= NBlist::NB_MAX"
	.text
	.p2align 4,,15
.globl gpuirr_set_list_
	.type	gpuirr_set_list_, @function
gpuirr_set_list_:
.LFB2402:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	movl	(%rsi), %eax
	movl	(%rdi), %r10d
	cmpl	$600, %eax
	jg	.L23
	decl	%r10d
	movslq	%r10d, %r10
	imulq	$2416, %r10, %r10
	addq	_ZL4list(%rip), %r10
	testl	%eax, %eax
	movl	%eax, 12(%r10)
	jle	.L9
	leaq	4(%rsi), %rdx
	leaq	16(%r10), %rdi
	movl	$4, %esi
	xorl	%ecx, %ecx
	vmovdqa	.LC2(%rip), %xmm0
	.p2align 4,,10
	.p2align 3
.L10:
	vmovdqu	(%rdx), %xmm2
	vmovdqu	16(%rdx), %xmm1
	vpsubd	%xmm0, %xmm2, %xmm2
	vpsubd	%xmm0, %xmm1, %xmm1
	vmovntdq	%xmm2, (%rdi)
	addl	$8, %ecx
	vmovntdq	%xmm1, 16(%r10,%rsi,4)
	addq	$32, %rdx
	addq	$32, %rdi
	addq	$8, %rsi
	cmpl	%ecx, %eax
	jg	.L10
.L9:
	leal	6(%rax), %r9d
	movl	%eax, %edx
	movl	_ZL4nmax(%rip), %esi
	decl	%edx
	cmovns	%edx, %r9d
	sarl	$3, %r9d
	leal	8(,%r9,8), %r9d
	cmpl	%r9d, %eax
	jge	.L17
	movl	%r9d, %r11d
	movslq	%eax, %rdx
	subl	%eax, %r11d
	leaq	4(%rdx), %rbp
	leaq	(%r10,%rbp,4), %r8
	andl	$15, %r8d
	shrq	$2, %r8
	negl	%r8d
	andl	$3, %r8d
	cmpl	%r11d, %r8d
	cmova	%r11d, %r8d
	testl	%r8d, %r8d
	je	.L12
	leaq	16(%r10,%rdx,4), %rcx
	movl	%eax, %edx
	.p2align 4,,10
	.p2align 3
.L13:
	incl	%edx
	movl	%esi, (%rcx)
	movl	%edx, %edi
	addq	$4, %rcx
	subl	%eax, %edi
	cmpl	%edi, %r8d
	ja	.L13
	cmpl	%r8d, %r11d
	je	.L17
	movl	%edx, %eax
.L12:
	subl	%r8d, %r11d
	movl	%r11d, %edi
	shrl	$2, %edi
	leal	0(,%rdi,4), %ebx
	testl	%ebx, %ebx
	je	.L14
	vmovd	%esi, %xmm1
	mov	%r8d, %r8d
	vpshufd	$0, %xmm1, %xmm0
	addq	%r8, %rbp
	xorl	%edx, %edx
	leaq	(%r10,%rbp,4), %rcx
	.p2align 4,,10
	.p2align 3
.L15:
	vmovdqa	%xmm0, (%rcx)
	incl	%edx
	addq	$16, %rcx
	cmpl	%edi, %edx
	jb	.L15
	addl	%ebx, %eax
	cmpl	%ebx, %r11d
	je	.L17
.L14:
	movslq	%eax, %rdx
	leaq	16(%r10,%rdx,4), %rdx
	.p2align 4,,10
	.p2align 3
.L16:
	movl	%esi, (%rdx)
	incl	%eax
	addq	$4, %rdx
	cmpl	%eax, %r9d
	jg	.L16
.L17:
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L23:
	.cfi_restore_state
	movl	$_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, %ecx
	movl	$433, %edx
	movl	$.LC0, %esi
	movl	$.LC1, %edi
	call	__assert_fail
	.cfi_endproc
.LFE2402:
	.size	gpuirr_set_list_, .-gpuirr_set_list_
	.p2align 4,,15
.globl gpuirr_set_jp_
	.type	gpuirr_set_jp_, @function
gpuirr_set_jp_:
.LFB2401:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	vmovsd	16(%rsi), %xmm1
	vmovsd	(%rsi), %xmm2
	vmovhpd	(%r9), %xmm1, %xmm0
	vmovsd	(%rdx), %xmm4
	vmovhpd	8(%rsi), %xmm2, %xmm1
	vmovhpd	8(%rdx), %xmm4, %xmm3
	vmovsd	(%rcx), %xmm5
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	vmovhpd	8(%rcx), %xmm5, %xmm4
	movq	16(%rbp), %rax
	vinsertf128	$0x1, %xmm0, %ymm1, %ymm1
	vmovsd	(%r8), %xmm6
	vcvtpd2psy	%ymm1, %xmm0
	vmovhpd	8(%r8), %xmm6, %xmm5
	vcvtps2pd	%xmm0, %ymm2
	vmovaps	%xmm0, %xmm0
	vsubpd	%ymm2, %ymm1, %ymm1
	vmovsd	16(%rdx), %xmm2
	vcvtpd2psy	%ymm1, %xmm1
	vinsertf128	$0x1, %xmm2, %ymm3, %ymm3
	vinsertf128	$0x1, %xmm1, %ymm0, %ymm0
	vmovsd	16(%rcx), %xmm2
	vmovddup	(%rax), %xmm1
	vinsertf128	$0x1, %xmm2, %ymm4, %ymm4
	movl	(%rdi), %eax
	vmovsd	16(%r8), %xmm2
	vcvtpd2psy	%ymm3, %xmm3
	vcvtpd2psy	%ymm4, %xmm4
	vinsertf128	$0x1, %xmm2, %ymm5, %ymm2
	vmovaps	%xmm3, %xmm3
	vcvtpd2psy	%ymm2, %xmm2
	vinsertf128	$0x1, %xmm4, %ymm3, %ymm3
	vmovaps	%xmm2, %xmm2
	decl	%eax
	vinsertf128	$0x1, %xmm1, %ymm2, %ymm2
	cltq
	leaq	(%rax,%rax,2), %rax
	salq	$5, %rax
	addq	_ZL4ptcl(%rip), %rax
	vmovaps	%ymm0, (%rax)
	vmovaps	%ymm3, 32(%rax)
	vmovaps	%ymm2, 64(%rax)
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE2401:
	.size	gpuirr_set_jp_, .-gpuirr_set_jp_
	.section	.rodata.str1.1
.LC3:
	.string	"gpuirr: it is already close\n"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC5:
	.string	"**************************** \n"
	.align 8
.LC6:
	.string	"Closing GPUIRR lib. CPU ver. \n"
	.section	.rodata.str1.1
.LC7:
	.string	"time grav  : %f sec\n"
.LC10:
	.string	"perf grav  : %f Gflops\n"
.LC11:
	.string	"perf grav  : %f usec\n"
.LC12:
	.string	"<#NB>      : %f \n"
	.text
	.p2align 4,,15
.globl gpuirr_close_
	.type	gpuirr_close_, @function
gpuirr_close_:
.LFB2400:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$136, %rsp
	.cfi_def_cfa_offset 144
	cmpb	$0, _ZL7is_open(%rip)
	je	.L35
	movq	_ZL4ptcl(%rip), %rdi
	call	free
	movq	$0, _ZL4ptcl(%rip)
	movq	_ZL4list(%rip), %rdi
	call	free
	movq	_ZL9num_inter(%rip), %rax
	movq	$0, _ZL4list(%rip)
	testq	%rax, %rax
	js	.L28
	vcvtsi2sdq	%rax, %xmm2, %xmm2
.L29:
	movq	_ZL9num_fcall(%rip), %rax
	vmovsd	_ZL9time_grav(%rip), %xmm1
	testq	%rax, %rax
	js	.L30
	vcvtsi2sdq	%rax, %xmm4, %xmm4
.L31:
	vdivsd	%xmm4, %xmm1, %xmm4
	movq	_ZL9num_steps(%rip), %rax
	vmulsd	.LC4(%rip), %xmm4, %xmm4
	testq	%rax, %rax
	js	.L32
	vcvtsi2sdq	%rax, %xmm3, %xmm3
.L33:
	vdivsd	%xmm3, %xmm2, %xmm3
	vmovsd	%xmm4, 32(%rsp)
	vmovsd	%xmm3, (%rsp)
	vmovsd	%xmm1, 64(%rsp)
	vmovsd	%xmm2, 96(%rsp)
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC5, %edi
	call	fwrite
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	vmovsd	_ZL9time_grav(%rip), %xmm0
	movl	$.LC7, %esi
	movq	stderr(%rip), %rdi
	movl	$1, %eax
	call	fprintf
	movq	stderr(%rip), %rsi
	movl	$10, %edi
	call	fputc
	vmovsd	96(%rsp), %xmm2
	vmovsd	64(%rsp), %xmm1
	vmulsd	.LC8(%rip), %xmm2, %xmm0
	movl	$.LC10, %esi
	vmulsd	.LC9(%rip), %xmm0, %xmm0
	movq	stderr(%rip), %rdi
	vdivsd	%xmm1, %xmm0, %xmm0
	movl	$1, %eax
	call	fprintf
	vmovsd	32(%rsp), %xmm4
	movl	$.LC11, %esi
	vmovapd	%xmm4, %xmm0
	movq	stderr(%rip), %rdi
	movl	$1, %eax
	call	fprintf
	vmovsd	(%rsp), %xmm3
	movl	$.LC12, %esi
	vmovapd	%xmm3, %xmm0
	movq	stderr(%rip), %rdi
	movl	$1, %eax
	call	fprintf
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC5, %edi
	call	fwrite
	movb	$0, _ZL7is_open(%rip)
	addq	$136, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L35:
	.cfi_restore_state
	movq	stderr(%rip), %rcx
	movl	$28, %edx
	movl	$1, %esi
	movl	$.LC3, %edi
	addq	$136, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	jmp	fwrite
	.p2align 4,,10
	.p2align 3
.L28:
	.cfi_restore_state
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm2, %xmm2
	vaddsd	%xmm2, %xmm2, %xmm2
	jmp	.L29
	.p2align 4,,10
	.p2align 3
.L32:
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm3, %xmm3
	vaddsd	%xmm3, %xmm3, %xmm3
	jmp	.L33
	.p2align 4,,10
	.p2align 3
.L30:
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm4, %xmm4
	vaddsd	%xmm4, %xmm4, %xmm4
	jmp	.L31
	.cfi_endproc
.LFE2400:
	.size	gpuirr_close_, .-gpuirr_close_
	.p2align 4,,15
	.type	_ZL11gpuirr_openii.omp_fn.1, @function
_ZL11gpuirr_openii.omp_fn.1:
.LFB2512:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	call	omp_get_num_threads
	movl	%eax, _ZL11num_threads(%rip)
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2512:
	.size	_ZL11gpuirr_openii.omp_fn.1, .-_ZL11gpuirr_openii.omp_fn.1
	.section	.rodata.str1.1
.LC13:
	.string	"gpuirr: it is already open\n"
.LC14:
	.string	"lmax <= 1 + NBlist::NB_MAX"
	.section	.rodata.str1.8
	.align 8
.LC15:
	.string	"Opening GPUIRR lib. AVX ver. \n"
	.section	.rodata.str1.1
.LC16:
	.string	" nmax = %d, lmax = %d\n"
	.section	.rodata.str1.8
	.align 8
.LC17:
	.string	"0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle))"
	.align 8
.LC19:
	.string	"0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist))"
	.text
	.p2align 4,,15
.globl gpuirr_open_
	.type	gpuirr_open_, @function
gpuirr_open_:
.LFB2399:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movl	(%rsi), %ebp
	movl	(%rdi), %ebx
	cmpb	$0, _ZL7is_open(%rip)
	jne	.L45
	cmpl	$601, %ebp
	jg	.L46
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC5, %edi
	leaq	8(%rsp), %r12
	call	fwrite
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC15, %edi
	call	fwrite
	movl	%ebp, %ecx
	movl	%ebx, %edx
	movl	$.LC16, %esi
	movq	stderr(%rip), %rdi
	xorl	%eax, %eax
	call	fprintf
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC5, %edi
	call	fwrite
	leal	1(%rbx), %eax
	movl	$64, %esi
	cltq
	movq	%r12, %rdi
	leaq	(%rax,%rax,2), %rbp
	salq	$5, %rbp
	movq	%rbp, %rdx
	call	posix_memalign
	testl	%eax, %eax
	jne	.L47
	movq	%rbp, %rdx
	movl	$255, %esi
	movq	8(%rsp), %rdi
	movslq	%ebx, %rbp
	call	memset
	movq	8(%rsp), %rdx
	leaq	0(%rbp,%rbp,2), %rax
	movq	%rdx, _ZL4ptcl(%rip)
	salq	$5, %rax
	vxorpd	%xmm0, %xmm0, %xmm0
	leaq	(%rdx,%rax), %rax
	imulq	$2416, %rbp, %rbp
	vmovapd	%xmm0, 80(%rax)
	movq	%rbp, %rdx
	vxorps	%xmm0, %xmm0, %xmm0
	movl	$64, %esi
	vmovaps	%xmm0, 16(%rax)
	vmovaps	%xmm0, 32(%rax)
	vmovaps	%xmm0, 48(%rax)
	vmovaps	%xmm0, 64(%rax)
	movq	%r12, %rdi
	vmovaps	.LC18(%rip), %xmm0
	vmovaps	%xmm0, (%rax)
	call	posix_memalign
	testl	%eax, %eax
	jne	.L48
	movq	%rbp, %rdx
	movq	8(%rsp), %rdi
	movl	$255, %esi
	call	memset
	movq	8(%rsp), %rax
	xorl	%edx, %edx
	movq	%rax, _ZL4list(%rip)
	xorl	%esi, %esi
	movl	%ebx, _ZL4nmax(%rip)
	movl	$_ZL11gpuirr_openii.omp_fn.1, %edi
	call	GOMP_parallel_start
	xorl	%edi, %edi
	call	_ZL11gpuirr_openii.omp_fn.1
	call	GOMP_parallel_end
	movq	$0, _ZL9time_grav(%rip)
	movq	$0, _ZL9num_steps(%rip)
	movq	$0, _ZL9num_fcall(%rip)
	movq	$0, _ZL9num_inter(%rip)
	movb	$1, _ZL7is_open(%rip)
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L45:
	.cfi_restore_state
	movq	stderr(%rip), %rcx
	movl	$27, %edx
	movl	$1, %esi
	movl	$.LC13, %edi
	call	fwrite
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
.L48:
	.cfi_restore_state
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$376, %edx
	movl	$.LC0, %esi
	movl	$.LC19, %edi
	call	__assert_fail
.L47:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$368, %edx
	movl	$.LC0, %esi
	movl	$.LC17, %edi
	call	__assert_fail
.L46:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$358, %edx
	movl	$.LC0, %esi
	movl	$.LC14, %edi
	call	__assert_fail
	.cfi_endproc
.LFE2399:
	.size	gpuirr_open_, .-gpuirr_open_
	.p2align 4,,15
	.type	_ZL11gpuirr_firriPdS_, @function
_ZL11gpuirr_firriPdS_:
.LFB2397:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movslq	%edi, %rdi
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	leaq	(%rdi,%rdi,2), %rcx
	andq	$-32, %rsp
	salq	$5, %rcx
	pushq	%rbx
	imulq	$2416, %rdi, %rbx
	.cfi_escape 0x10,0x3,0x7,0x76,0x0,0x9,0xe0,0x1a,0x38,0x1c
	subq	$480, %rsp
	vmovapd	_ZL8vec_tnow(%rip), %xmm0
	addq	_ZL4list(%rip), %rbx
	vmovapd	%xmm0, 440(%rsp)
	movq	_ZL4ptcl(%rip), %rax
	vmovapd	440(%rsp), %xmm1
	leaq	(%rax,%rcx), %rcx
	movl	12(%rbx), %r11d
	vmovapd	80(%rcx), %xmm0
	testl	%r11d, %r11d
	vmovddup	%xmm0, %xmm0
	vsubpd	%xmm0, %xmm1, %xmm4
	vcvtpd2psx	%xmm4, %xmm4
	vshufps	$0, %xmm4, %xmm4, %xmm0
	vshufps	$85, %xmm4, %xmm4, %xmm4
	vmovaps	%xmm0, %xmm0
	vinsertf128	$0x1, %xmm4, %ymm0, %ymm4
	vmovaps	16(%rcx), %xmm0
	vaddps	%ymm4, %ymm4, %ymm7
	vmovaps	%xmm0, %xmm2
	vinsertf128	$0x1, %xmm0, %ymm2, %ymm2
	vmovaps	32(%rcx), %xmm0
	vmovaps	%xmm0, %xmm1
	vinsertf128	$0x1, %xmm0, %ymm1, %ymm0
	vmovaps	48(%rcx), %xmm1
	vmovaps	%xmm1, %xmm5
	vinsertf128	$0x1, %xmm1, %ymm5, %ymm5
	vmovaps	64(%rcx), %xmm1
	vmovaps	%xmm1, %xmm6
	vinsertf128	$0x1, %xmm1, %ymm6, %ymm6
	vmovaps	(%rcx), %xmm1
	vmovaps	%xmm1, %xmm3
	vinsertf128	$0x1, %xmm1, %ymm3, %ymm3
	vmulps	%ymm6, %ymm4, %ymm1
	vaddps	%ymm1, %ymm5, %ymm1
	vmulps	%ymm1, %ymm4, %ymm1
	vaddps	%ymm1, %ymm0, %ymm1
	vmulps	%ymm1, %ymm4, %ymm1
	vaddps	%ymm1, %ymm2, %ymm2
	vmovaps	.LC21(%rip), %ymm1
	vmulps	%ymm1, %ymm4, %ymm4
	vmulps	%ymm6, %ymm4, %ymm6
	vaddps	%ymm6, %ymm5, %ymm5
	vshufps	$0, %ymm2, %ymm2, %ymm6
	vmulps	%ymm5, %ymm7, %ymm4
	vmovaps	%ymm6, 72(%rsp)
	vaddps	%ymm4, %ymm0, %ymm0
	vshufps	$85, %ymm3, %ymm3, %ymm5
	vshufps	$0, %ymm3, %ymm3, %ymm4
	vshufps	$85, %ymm2, %ymm2, %ymm7
	vshufps	$170, %ymm3, %ymm3, %ymm3
	vshufps	$170, %ymm2, %ymm2, %ymm2
	vmovaps	%ymm3, 104(%rsp)
	vmovaps	%ymm2, 8(%rsp)
	vshufps	$85, %ymm0, %ymm0, %ymm3
	vshufps	$0, %ymm0, %ymm0, %ymm2
	vmovaps	%ymm4, 168(%rsp)
	vshufps	$170, %ymm0, %ymm0, %ymm0
	vmovaps	%ymm5, 136(%rsp)
	vmovaps	%ymm7, 40(%rsp)
	vmovaps	%ymm2, -24(%rsp)
	vmovaps	%ymm3, -56(%rsp)
	vmovaps	%ymm0, -88(%rsp)
	jle	.L57
	movq	%rbx, %r10
	movq	%rbx, %r8
	xorl	%edi, %edi
	.p2align 4,,10
	.p2align 3
.L52:
	movslq	16(%r8), %rcx
	incl	%edi
	leaq	(%rcx,%rcx,2), %rcx
	addq	$4, %r8
	salq	$5, %rcx
	cmpl	%edi, %r11d
	leaq	(%rax,%rcx), %rcx
	prefetcht0	(%rcx)
	prefetcht0	64(%rcx)
	jg	.L52
	vxorps	%xmm5, %xmm5, %xmm5
	xorl	%r9d, %r9d
	vmovaps	%ymm5, 360(%rsp)
	vmovaps	%ymm5, 200(%rsp)
	vmovaps	%ymm5, 232(%rsp)
	vmovaps	%ymm5, 328(%rsp)
	vmovaps	%ymm5, 296(%rsp)
	vmovaps	%ymm5, 264(%rsp)
	.p2align 4,,10
	.p2align 3
.L53:
	movslq	16(%r10), %r8
	vmovapd	440(%rsp), %xmm6
	vmovapd	440(%rsp), %xmm7
	movslq	%r9d, %rcx
	leaq	(%r8,%r8,2), %r8
	leaq	16(%rbx,%rcx,4), %rcx
	salq	$5, %r8
	movslq	4(%rcx), %rdi
	leaq	(%rax,%r8), %r8
	leaq	(%rdi,%rdi,2), %rdi
	vmovapd	80(%r8), %xmm0
	vmovaps	16(%r8), %xmm8
	vmovaps	48(%r8), %xmm5
	vmovaps	(%r8), %xmm12
	salq	$5, %rdi
	vmovaps	%xmm8, %xmm8
	leaq	(%rax,%rdi), %rdi
	vmovaps	%xmm5, %xmm5
	vmovhpd	80(%rdi), %xmm0, %xmm0
	vinsertf128	$0x1, 48(%rdi), %ymm5, %ymm5
	vinsertf128	$0x1, 16(%rdi), %ymm8, %ymm8
	vsubpd	%xmm0, %xmm6, %xmm0
	vmovaps	%xmm12, %xmm12
	vmovaps	64(%r8), %xmm6
	vinsertf128	$0x1, (%rdi), %ymm12, %ymm12
	vcvtpd2psx	%xmm0, %xmm0
	vmovaps	%xmm6, %xmm6
	vshufps	$0, %xmm0, %xmm0, %xmm2
	vinsertf128	$0x1, 64(%rdi), %ymm6, %ymm6
	vmovaps	%xmm2, %xmm2
	vshufps	$85, %xmm0, %xmm0, %xmm0
	addl	$8, %r9d
	vinsertf128	$0x1, %xmm0, %ymm2, %ymm0
	addq	$32, %r10
	vmovaps	32(%r8), %xmm2
	vmulps	%ymm6, %ymm0, %ymm3
	movslq	8(%rcx), %r8
	vaddps	%ymm3, %ymm5, %ymm3
	vmovaps	%xmm2, %xmm2
	vmulps	%ymm3, %ymm0, %ymm3
	vinsertf128	$0x1, 32(%rdi), %ymm2, %ymm2
	leaq	(%r8,%r8,2), %r8
	movslq	12(%rcx), %rdi
	vaddps	%ymm3, %ymm2, %ymm3
	salq	$5, %r8
	vmulps	%ymm3, %ymm0, %ymm3
	leaq	(%rax,%r8), %r8
	vaddps	%ymm3, %ymm8, %ymm8
	vmovaps	64(%r8), %xmm11
	vaddps	%ymm0, %ymm0, %ymm3
	vmovaps	(%r8), %xmm10
	vmulps	%ymm1, %ymm0, %ymm0
	leaq	(%rdi,%rdi,2), %rdi
	vmulps	%ymm6, %ymm0, %ymm0
	salq	$5, %rdi
	vmovaps	16(%r8), %xmm6
	leaq	(%rax,%rdi), %rdi
	vaddps	%ymm0, %ymm5, %ymm0
	vmovaps	%xmm6, %xmm6
	vmulps	%ymm0, %ymm3, %ymm0
	vinsertf128	$0x1, 16(%rdi), %ymm6, %ymm6
	vaddps	%ymm0, %ymm2, %ymm4
	vmovaps	%xmm11, %xmm11
	vmovapd	80(%r8), %xmm0
	vinsertf128	$0x1, 64(%rdi), %ymm11, %ymm11
	vmovhpd	80(%rdi), %xmm0, %xmm0
	vmovaps	%xmm10, %xmm10
	vsubpd	%xmm0, %xmm7, %xmm0
	vinsertf128	$0x1, (%rdi), %ymm10, %ymm10
	vcvtpd2psx	%xmm0, %xmm0
	vmovaps	48(%r8), %xmm7
	vshufps	$85, %xmm0, %xmm0, %xmm3
	vshufps	$0, %xmm0, %xmm0, %xmm2
	vmovaps	%xmm7, %xmm7
	vmovaps	%xmm2, %xmm0
	vinsertf128	$0x1, 48(%rdi), %ymm7, %ymm7
	vinsertf128	$0x1, %xmm3, %ymm0, %ymm0
	vmovaps	32(%r8), %xmm3
	vmulps	%ymm11, %ymm0, %ymm2
	vaddps	%ymm0, %ymm0, %ymm5
	movslq	16(%rcx), %r8
	vaddps	%ymm2, %ymm7, %ymm2
	vmovaps	%xmm3, %xmm3
	vmulps	%ymm2, %ymm0, %ymm2
	vinsertf128	$0x1, 32(%rdi), %ymm3, %ymm3
	leaq	(%r8,%r8,2), %r8
	movslq	20(%rcx), %rdi
	vaddps	%ymm2, %ymm3, %ymm2
	salq	$5, %r8
	vmulps	%ymm2, %ymm0, %ymm2
	leaq	(%rax,%r8), %r8
	vaddps	%ymm2, %ymm6, %ymm6
	vmovaps	48(%r8), %xmm9
	vmovaps	64(%r8), %xmm13
	vmulps	%ymm1, %ymm0, %ymm0
	leaq	(%rdi,%rdi,2), %rdi
	vmulps	%ymm11, %ymm0, %ymm0
	salq	$5, %rdi
	vaddps	%ymm0, %ymm7, %ymm0
	leaq	(%rax,%rdi), %rdi
	vmulps	%ymm0, %ymm5, %ymm0
	vmovaps	16(%r8), %xmm7
	vaddps	%ymm0, %ymm3, %ymm2
	vmovaps	(%r8), %xmm11
	vmovapd	440(%rsp), %xmm3
	vmovapd	80(%r8), %xmm0
	vmovaps	%xmm7, %xmm7
	vmovhpd	80(%rdi), %xmm0, %xmm0
	vinsertf128	$0x1, 16(%rdi), %ymm7, %ymm7
	vsubpd	%xmm0, %xmm3, %xmm0
	vmovaps	%xmm9, %xmm9
	vcvtpd2psx	%xmm0, %xmm0
	vinsertf128	$0x1, 48(%rdi), %ymm9, %ymm9
	vshufps	$85, %xmm0, %xmm0, %xmm5
	vshufps	$0, %xmm0, %xmm0, %xmm3
	vmovaps	%xmm13, %xmm13
	vmovaps	%xmm3, %xmm0
	vinsertf128	$0x1, 64(%rdi), %ymm13, %ymm13
	vinsertf128	$0x1, %xmm5, %ymm0, %ymm0
	vmovaps	%xmm11, %xmm11
	vmulps	%ymm13, %ymm0, %ymm3
	vinsertf128	$0x1, (%rdi), %ymm11, %ymm11
	vaddps	%ymm0, %ymm0, %ymm14
	vmovaps	32(%r8), %xmm5
	vaddps	%ymm3, %ymm9, %ymm3
	vmovaps	%xmm5, %xmm5
	vmulps	%ymm3, %ymm0, %ymm3
	vinsertf128	$0x1, 32(%rdi), %ymm5, %ymm5
	vaddps	%ymm3, %ymm5, %ymm3
	vmulps	%ymm3, %ymm0, %ymm3
	vmulps	%ymm1, %ymm0, %ymm0
	vaddps	%ymm3, %ymm7, %ymm7
	vmulps	%ymm13, %ymm0, %ymm0
	vaddps	%ymm0, %ymm9, %ymm0
	vmulps	%ymm0, %ymm14, %ymm0
	vaddps	%ymm0, %ymm5, %ymm0
	vmovaps	%ymm0, -120(%rsp)
	movslq	28(%rcx), %rdi
	vmovapd	440(%rsp), %xmm5
	movslq	24(%rcx), %rcx
	leaq	(%rdi,%rdi,2), %rdi
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rdi
	salq	$5, %rcx
	leaq	(%rax,%rdi), %rdi
	leaq	(%rax,%rcx), %rcx
	cmpl	%r9d, %r11d
	vmovapd	80(%rcx), %xmm0
	vmovaps	32(%rcx), %xmm13
	vmovhpd	80(%rdi), %xmm0, %xmm0
	vmovaps	48(%rcx), %xmm14
	vsubpd	%xmm0, %xmm5, %xmm0
	vmovaps	64(%rcx), %xmm15
	vcvtpd2psx	%xmm0, %xmm0
	vmovaps	%xmm13, %xmm13
	vshufps	$85, %xmm0, %xmm0, %xmm9
	vinsertf128	$0x1, 32(%rdi), %ymm13, %ymm13
	vshufps	$0, %xmm0, %xmm0, %xmm5
	vmovaps	%xmm14, %xmm14
	vmovaps	%xmm5, %xmm0
	vinsertf128	$0x1, 48(%rdi), %ymm14, %ymm14
	vinsertf128	$0x1, %xmm9, %ymm0, %ymm0
	vmovaps	16(%rcx), %xmm5
	vmovaps	(%rcx), %xmm9
	vmovaps	%xmm5, %xmm5
	vmovaps	%xmm15, %xmm15
	vinsertf128	$0x1, 16(%rdi), %ymm5, %ymm3
	vinsertf128	$0x1, 64(%rdi), %ymm15, %ymm15
	vmovaps	%xmm9, %xmm9
	vmulps	%ymm15, %ymm0, %ymm5
	vinsertf128	$0x1, (%rdi), %ymm9, %ymm9
	vaddps	%ymm5, %ymm14, %ymm5
	vmulps	%ymm5, %ymm0, %ymm5
	vaddps	%ymm5, %ymm13, %ymm5
	vmulps	%ymm5, %ymm0, %ymm5
	vaddps	%ymm3, %ymm5, %ymm5
	vaddps	%ymm0, %ymm0, %ymm3
	vmulps	%ymm1, %ymm0, %ymm0
	vmulps	%ymm15, %ymm0, %ymm0
	vaddps	%ymm0, %ymm14, %ymm0
	vunpcklps	%ymm9, %ymm10, %ymm14
	vmulps	%ymm3, %ymm0, %ymm0
	vunpckhps	%ymm9, %ymm10, %ymm9
	vaddps	%ymm0, %ymm13, %ymm13
	vunpckhps	-120(%rsp), %ymm4, %ymm3
	vunpcklps	%ymm11, %ymm12, %ymm0
	vunpckhps	%ymm11, %ymm12, %ymm11
	vunpcklps	%ymm14, %ymm0, %ymm15
	vunpcklps	%ymm9, %ymm11, %ymm10
	vunpckhps	%ymm14, %ymm0, %ymm14
	vunpckhps	%ymm9, %ymm11, %ymm11
	vunpcklps	%ymm7, %ymm8, %ymm0
	vunpcklps	%ymm5, %ymm6, %ymm9
	vsubps	168(%rsp), %ymm15, %ymm15
	vunpckhps	%ymm5, %ymm6, %ymm5
	vunpcklps	%ymm9, %ymm0, %ymm12
	vsubps	136(%rsp), %ymm14, %ymm14
	vunpckhps	%ymm9, %ymm0, %ymm9
	vsubps	72(%rsp), %ymm12, %ymm12
	vunpcklps	-120(%rsp), %ymm4, %ymm0
	vaddps	%ymm12, %ymm15, %ymm15
	vsubps	40(%rsp), %ymm9, %ymm9
	vsubps	104(%rsp), %ymm10, %ymm10
	vaddps	%ymm9, %ymm14, %ymm14
	vunpckhps	%ymm7, %ymm8, %ymm7
	vunpcklps	%ymm5, %ymm7, %ymm7
	vunpcklps	%ymm13, %ymm2, %ymm5
	vsubps	8(%rsp), %ymm7, %ymm7
	vunpcklps	%ymm5, %ymm0, %ymm6
	vaddps	%ymm7, %ymm10, %ymm7
	vsubps	-24(%rsp), %ymm6, %ymm6
	vunpckhps	%ymm13, %ymm2, %ymm13
	vunpckhps	%ymm5, %ymm0, %ymm5
	vmulps	%ymm15, %ymm15, %ymm2
	vsubps	-56(%rsp), %ymm5, %ymm5
	vmulps	%ymm14, %ymm14, %ymm0
	vunpcklps	%ymm13, %ymm3, %ymm3
	vaddps	%ymm0, %ymm2, %ymm2
	vsubps	-88(%rsp), %ymm3, %ymm3
	vmulps	%ymm7, %ymm7, %ymm0
	vaddps	%ymm0, %ymm2, %ymm2
	vrsqrtps	%ymm2, %ymm0
	vmulps	.LC22(%rip), %ymm0, %ymm4
	vmulps	%ymm0, %ymm2, %ymm2
	vmulps	%ymm2, %ymm0, %ymm2
	vaddps	.LC23(%rip), %ymm2, %ymm9
	vmulps	%ymm5, %ymm14, %ymm2
	vmulps	%ymm9, %ymm4, %ymm9
	vmulps	%ymm6, %ymm15, %ymm4
	vmulps	%ymm9, %ymm9, %ymm0
	vaddps	%ymm2, %ymm4, %ymm4
	vmulps	%ymm11, %ymm9, %ymm11
	vmulps	%ymm3, %ymm7, %ymm2
	vaddps	%ymm2, %ymm4, %ymm4
	vmulps	.LC23(%rip), %ymm0, %ymm2
	vmulps	%ymm11, %ymm0, %ymm0
	vmulps	%ymm2, %ymm4, %ymm2
	vmulps	%ymm0, %ymm15, %ymm4
	vmulps	%ymm2, %ymm15, %ymm15
	vaddps	296(%rsp), %ymm4, %ymm4
	vaddps	%ymm15, %ymm6, %ymm6
	vmovaps	%ymm4, 296(%rsp)
	vmulps	%ymm6, %ymm0, %ymm6
	vmulps	%ymm0, %ymm14, %ymm4
	vmulps	%ymm2, %ymm14, %ymm14
	vaddps	328(%rsp), %ymm4, %ymm4
	vaddps	%ymm14, %ymm5, %ymm5
	vmovaps	%ymm4, 328(%rsp)
	vmulps	%ymm5, %ymm0, %ymm5
	vmulps	%ymm0, %ymm7, %ymm4
	vmulps	%ymm2, %ymm7, %ymm2
	vaddps	360(%rsp), %ymm4, %ymm4
	vaddps	%ymm2, %ymm3, %ymm3
	vmovaps	%ymm4, 360(%rsp)
	vmulps	%ymm3, %ymm0, %ymm3
	vaddps	200(%rsp), %ymm6, %ymm6
	vaddps	232(%rsp), %ymm5, %ymm5
	vaddps	264(%rsp), %ymm3, %ymm3
	vmovaps	%ymm6, 200(%rsp)
	vmovaps	%ymm5, 232(%rsp)
	vmovaps	%ymm3, 264(%rsp)
	jg	.L53
.L51:
	vmovaps	296(%rsp), %ymm4
	vmovaps	328(%rsp), %ymm5
	vmovaps	360(%rsp), %ymm6
	vmovaps	200(%rsp), %ymm7
	vextractf128	$0x0, %ymm4, %xmm1
	vextractf128	$0x1, %ymm4, %xmm0
	vcvtps2pd	%xmm1, %ymm1
	vextractf128	$0x0, %ymm5, %xmm2
	vcvtps2pd	%xmm0, %ymm0
	vcvtps2pd	%xmm2, %ymm2
	vaddpd	%ymm0, %ymm1, %ymm0
	vhaddpd	%ymm0, %ymm0, %ymm0
	vextractf128	$0x0, %ymm0, %xmm1
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm1, %xmm1
	vextractf128	$0x1, %ymm5, %xmm0
	vmovlpd	%xmm1, (%rsi)
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm2, %ymm0
	vextractf128	$0x0, %ymm6, %xmm2
	vhaddpd	%ymm0, %ymm0, %ymm0
	vcvtps2pd	%xmm2, %ymm2
	vextractf128	$0x0, %ymm0, %xmm3
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm3, %xmm3
	vextractf128	$0x1, %ymm6, %xmm0
	vmovlpd	%xmm3, 8(%rsi)
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm2, %ymm0
	vextractf128	$0x0, %ymm7, %xmm2
	vhaddpd	%ymm0, %ymm0, %ymm0
	vcvtps2pd	%xmm2, %ymm2
	vextractf128	$0x0, %ymm0, %xmm5
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm5, %xmm5
	vextractf128	$0x1, %ymm7, %xmm0
	vmovlpd	%xmm5, 16(%rsi)
	vmovaps	264(%rsp), %ymm7
	vcvtps2pd	%xmm0, %ymm0
	vextractf128	$0x0, %ymm7, %xmm6
	vaddpd	%ymm0, %ymm2, %ymm0
	vcvtps2pd	%xmm6, %ymm6
	vhaddpd	%ymm0, %ymm0, %ymm0
	vextractf128	$0x0, %ymm0, %xmm4
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm4, %xmm4
	vmovaps	232(%rsp), %ymm0
	vmovlpd	%xmm4, (%rdx)
	vextractf128	$0x0, %ymm0, %xmm2
	vextractf128	$0x1, %ymm0, %xmm0
	vcvtps2pd	%xmm2, %ymm2
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm2, %ymm0
	vhaddpd	%ymm0, %ymm0, %ymm0
	vextractf128	$0x0, %ymm0, %xmm2
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm2, %xmm2
	vextractf128	$0x1, %ymm7, %xmm0
	vmovlpd	%xmm2, 8(%rdx)
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm6, %ymm0
	vhaddpd	%ymm0, %ymm0, %ymm0
	vextractf128	$0x0, %ymm0, %xmm6
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm6, %xmm0
	vmovlpd	%xmm0, 16(%rdx)
	addq	$480, %rsp
	popq	%rbx
	leave
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L57:
	.cfi_restore_state
	vxorps	%xmm4, %xmm4, %xmm4
	vmovaps	%ymm4, 360(%rsp)
	vmovaps	%ymm4, 200(%rsp)
	vmovaps	%ymm4, 232(%rsp)
	vmovaps	%ymm4, 328(%rsp)
	vmovaps	%ymm4, 296(%rsp)
	vmovaps	%ymm4, 264(%rsp)
	jmp	.L51
	.cfi_endproc
.LFE2397:
	.size	_ZL11gpuirr_firriPdS_, .-_ZL11gpuirr_firriPdS_
	.p2align 4,,15
	.type	_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0, @function
_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0:
.LFB2511:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movl	$1, %ecx
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movl	$1, %edx
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	xorl	%r14d, %r14d
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	xorl	%r13d, %r13d
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, %rbx
	subq	$40, %rsp
	.cfi_def_cfa_offset 96
	movslq	24(%rdi), %rsi
	leaq	16(%rsp), %rbp
	xorl	%edi, %edi
	leaq	24(%rsp), %r12
	movq	%rbp, %r9
	movq	%r12, %r8
	call	GOMP_loop_guided_start
	testb	%al, %al
	je	.L60
	.p2align 4,,10
	.p2align 3
.L64:
	movl	24(%rsp), %eax
	movl	16(%rsp), %edx
	movslq	%eax, %r14
	movl	%edx, 12(%rsp)
	leaq	(%r14,%r14,2), %r15
	movq	(%rbx), %rcx
	salq	$3, %r15
	salq	$2, %r14
	.p2align 4,,10
	.p2align 3
.L61:
	movl	(%rcx,%r14), %edi
	movq	%r15, %rdx
	movq	%r15, %rsi
	addq	16(%rbx), %rdx
	addq	8(%rbx), %rsi
	decl	%edi
	movl	%eax, (%rsp)
	addq	$24, %r15
	call	_ZL11gpuirr_firriPdS_
	movq	(%rbx), %rcx
	movq	_ZL4list(%rip), %rdx
	movslq	(%rcx,%r14), %rsi
	movl	(%rsp), %eax
	decq	%rsi
	incl	%eax
	imulq	$2416, %rsi, %rsi
	addq	$4, %r14
	addl	12(%rsi,%rdx), %r13d
	cmpl	%eax, 12(%rsp)
	jg	.L61
	movq	%rbp, %rsi
	movq	%r12, %rdi
	call	GOMP_loop_guided_next
	testb	%al, %al
	jne	.L64
	movl	%r13d, %r14d
.L60:
	call	GOMP_loop_end_nowait
	lock addl	%r14d, 28(%rbx)
	addq	$40, %rsp
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
.LFE2511:
	.size	_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0, .-_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0
	.p2align 4,,15
.globl gpuirr_firr_vec_
	.type	gpuirr_firr_vec_, @function
gpuirr_firr_vec_:
.LFB2405:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	movq	%rbx, -40(%rsp)
	movq	%r12, -24(%rsp)
	movq	%r13, -16(%rsp)
	movq	%r14, -8(%rsp)
	movq	%rdx, %r13
	.cfi_offset 14, -16
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	.cfi_offset 3, -48
	movq	%rcx, %r14
	movq	%rbp, -32(%rsp)
	movq	%rsi, %rbx
	subq	$136, %rsp
	.cfi_def_cfa_offset 144
	.cfi_offset 6, -40
	xorl	%esi, %esi
	movslq	(%rdi), %rbp
	leaq	80(%rsp), %r12
	movq	%r12, %rdi
	call	gettimeofday
	vcvtsi2sdq	80(%rsp), %xmm2, %xmm2
	vmovsd	.LC24(%rip), %xmm0
	xorl	%edx, %edx
	vmovsd	%xmm0, (%rsp)
	vcvtsi2sdq	88(%rsp), %xmm1, %xmm1
	movl	%ebp, 72(%rsp)
	vmulsd	%xmm0, %xmm1, %xmm1
	movq	%rbx, 48(%rsp)
	vaddsd	%xmm1, %xmm2, %xmm1
	leaq	48(%rsp), %rbx
	vmovsd	%xmm1, 40(%rsp)
	movq	%rbx, %rsi
	movq	%r13, 56(%rsp)
	movq	%r14, 64(%rsp)
	movl	$0, 76(%rsp)
	movl	$_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0, %edi
	call	GOMP_parallel_start
	movq	%rbx, %rdi
	call	_ZL15gpuirr_firr_veciPKiPA3_dS2_.omp_fn.0
	call	GOMP_parallel_end
	movslq	76(%rsp), %rax
	movq	%r12, %rdi
	addq	%rax, _ZL9num_inter(%rip)
	xorl	%esi, %esi
	call	gettimeofday
	vcvtsi2sdq	80(%rsp), %xmm1, %xmm1
	vcvtsi2sdq	88(%rsp), %xmm2, %xmm2
	vmovsd	(%rsp), %xmm0
	addq	%rbp, _ZL9num_steps(%rip)
	vmulsd	%xmm0, %xmm2, %xmm0
	incq	_ZL9num_fcall(%rip)
	vaddsd	%xmm0, %xmm1, %xmm0
	movq	96(%rsp), %rbx
	vsubsd	40(%rsp), %xmm0, %xmm0
	movq	104(%rsp), %rbp
	vaddsd	_ZL9time_grav(%rip), %xmm0, %xmm0
	movq	112(%rsp), %r12
	vmovsd	%xmm0, _ZL9time_grav(%rip)
	movq	120(%rsp), %r13
	movq	128(%rsp), %r14
	addq	$136, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2405:
	.size	gpuirr_firr_vec_, .-gpuirr_firr_vec_
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.local	_ZL4list
	.comm	_ZL4list,8,8
	.local	_ZL9num_inter
	.comm	_ZL9num_inter,8,8
	.local	_ZL9time_grav
	.comm	_ZL9time_grav,8,8
	.local	_ZL9num_fcall
	.comm	_ZL9num_fcall,8,8
	.local	_ZL9num_steps
	.comm	_ZL9num_steps,8,8
	.local	_ZL4ptcl
	.comm	_ZL4ptcl,8,8
	.local	_ZL8vec_tnow
	.comm	_ZL8vec_tnow,16,16
	.local	_ZL4nmax
	.comm	_ZL4nmax,4,4
	.section	.rodata
	.align 32
	.type	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, 43
_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_set_list(int, int, const int*)"
	.local	_ZL7is_open
	.comm	_ZL7is_open,1,1
	.local	_ZL11num_threads
	.comm	_ZL11num_threads,4,4
	.align 16
	.type	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, 27
_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_open(int, int)"
	.weakref	_ZL20__gthrw_pthread_oncePiPFvvE,pthread_once
	.weakref	_ZL27__gthrw_pthread_getspecificj,pthread_getspecific
	.weakref	_ZL27__gthrw_pthread_setspecificjPKv,pthread_setspecific
	.weakref	_ZL22__gthrw_pthread_createPmPK14pthread_attr_tPFPvS3_ES3_,pthread_create
	.weakref	_ZL20__gthrw_pthread_joinmPPv,pthread_join
	.weakref	_ZL21__gthrw_pthread_equalmm,pthread_equal
	.weakref	_ZL20__gthrw_pthread_selfv,pthread_self
	.weakref	_ZL22__gthrw_pthread_detachm,pthread_detach
	.weakref	_ZL22__gthrw_pthread_cancelm,pthread_cancel
	.weakref	_ZL19__gthrw_sched_yieldv,sched_yield
	.weakref	_ZL26__gthrw_pthread_mutex_lockP15pthread_mutex_t,pthread_mutex_lock
	.weakref	_ZL29__gthrw_pthread_mutex_trylockP15pthread_mutex_t,pthread_mutex_trylock
	.weakref	_ZL31__gthrw_pthread_mutex_timedlockP15pthread_mutex_tPK8timespec,pthread_mutex_timedlock
	.weakref	_ZL28__gthrw_pthread_mutex_unlockP15pthread_mutex_t,pthread_mutex_unlock
	.weakref	_ZL26__gthrw_pthread_mutex_initP15pthread_mutex_tPK19pthread_mutexattr_t,pthread_mutex_init
	.weakref	_ZL29__gthrw_pthread_mutex_destroyP15pthread_mutex_t,pthread_mutex_destroy
	.weakref	_ZL30__gthrw_pthread_cond_broadcastP14pthread_cond_t,pthread_cond_broadcast
	.weakref	_ZL27__gthrw_pthread_cond_signalP14pthread_cond_t,pthread_cond_signal
	.weakref	_ZL25__gthrw_pthread_cond_waitP14pthread_cond_tP15pthread_mutex_t,pthread_cond_wait
	.weakref	_ZL30__gthrw_pthread_cond_timedwaitP14pthread_cond_tP15pthread_mutex_tPK8timespec,pthread_cond_timedwait
	.weakref	_ZL28__gthrw_pthread_cond_destroyP14pthread_cond_t,pthread_cond_destroy
	.weakref	_ZL26__gthrw_pthread_key_createPjPFvPvE,pthread_key_create
	.weakref	_ZL26__gthrw_pthread_key_deletej,pthread_key_delete
	.weakref	_ZL30__gthrw_pthread_mutexattr_initP19pthread_mutexattr_t,pthread_mutexattr_init
	.weakref	_ZL33__gthrw_pthread_mutexattr_settypeP19pthread_mutexattr_ti,pthread_mutexattr_settype
	.weakref	_ZL33__gthrw_pthread_mutexattr_destroyP19pthread_mutexattr_t,pthread_mutexattr_destroy
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC2:
	.long	1
	.long	1
	.long	1
	.long	1
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC4:
	.long	0
	.long	1093567616
	.align 8
.LC8:
	.long	0
	.long	1078853632
	.align 8
.LC9:
	.long	3894859413
	.long	1041313291
	.section	.rodata.cst16
	.align 16
.LC18:
	.long	1132396544
	.long	1132396544
	.long	1132396544
	.long	0
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC21:
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.align 32
.LC22:
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.align 32
.LC23:
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.section	.rodata.cst8
	.align 8
.LC24:
	.long	2696277389
	.long	1051772663
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-3)"
	.section	.note.GNU-stack,"",@progbits
