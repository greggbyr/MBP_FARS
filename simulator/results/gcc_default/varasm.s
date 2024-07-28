gcc2_compiled.:
___gnu_compiled_c:
.data
	.align 4
_in_section:
	.word	0
.text
	.align 8
LC0:
	.ascii "%s\12\0"
	.align 8
LC1:
	.ascii ".text\0"
	.align 4
	.global _text_section
	.proc	020
_text_section:
	!#PROLOGUE# 0
	save %sp,-112,%sp
	!#PROLOGUE# 1
	sethi %hi(_in_section),%l0
	ld [%l0+%lo(_in_section)],%o0
	cmp %o0,1
	be L4
	sethi %hi(_asm_out_file),%o0
	ld [%o0+%lo(_asm_out_file)],%o0
	sethi %hi(LC0),%o1
	or %o1,%lo(LC0),%o1
	sethi %hi(LC1),%o2
	call _fprintf,0
	or %o2,%lo(LC1),%o2
	mov 1,%o0
	st %o0,[%l0+%lo(_in_section)]
L4:
	ret
	restore
	.align 8
LC2:
	.ascii ".data\0"
	.align 4
	.global _data_section
	.proc	020
_data_section:
	!#PROLOGUE# 0
	save %sp,-112,%sp
	!#PROLOGUE# 1
	sethi %hi(_in_section),%o0
	ld [%o0+%lo(_in_section)],%o0
	cmp %o0,2
	be L10
	sethi %hi(_asm_out_file),%o0
	ld [%o0+%lo(_asm_out_file)],%o0
	sethi %hi(LC0),%o1
	or %o1,%lo(LC0),%o1
	sethi %hi(LC2),%o2
	call _fprintf,0
	or %o2,%lo(LC2),%o2
	mov 2,%o0
	sethi %hi(_in_section),%o1
	st %o0,[%o1+%lo(_in_section)]
L10:
	ret
	restore
	.align 4
	.global _make_function_rtl
	.proc	020
_make_function_rtl:
	!#PROLOGUE# 0
	save %sp,-112,%sp
	!#PROLOGUE# 1
	ld [%i0+64],%o0
	cmp %o0,0
	bne,a L17
	mov 1,%o0
	ld [%i0+36],%o1
	ld [%o1+20],%o2
	mov 39,%o0
	call _gen_rtx,0
	mov 4,%o1
	mov %o0,%o2
	ld [%i0+28],%o1
	call _gen_rtx,0
	mov 37,%o0
	st %o0,[%i0+64]
	mov 1,%o0
L17:
	sethi %hi(_function_defined),%o1
	st %o0,[%o1+%lo(_function_defined)]
	ret
	restore
