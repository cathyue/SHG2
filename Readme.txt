Major programs: 

A->B, means A depends on B, A is the major program, or caller


dispersion calculation: detuning->ome_lp

detuning: given l, p, n, R, calculate (w2-w)/(w2) and the off-resonance factor (see test.pdf), with figure;

ome_lp->n_lam
ome_lp: given l,p,n, R, calculate omega;

n_lam: calculate n_lam(microns), Kozyreff Appendix D

-------------------------------

intrinsic parameters defination and calculation: para->(sample from detuning), Aij_cal, Int1_cal, Gm2, jl, hl

para: define some intrinsic constants, calculate Bij and kmm

Aij_cal->Int1_cal, phil, kail
Aij_cal: given R, k1, n1, l1, k2, n2, l2, calculate Aij with 1d integral.
my notebook, calculate A11, A12...

Int1_cal->phil, kail
Int1_cal: [ Int1 ] = Int1_cal( R, k, n, l ), 
my notebook, ->only taking into account m=l>>1

Gm2->jl, hl
modified Gm, in my notebook.

jl, hl, phil, kail: spherical/Ricatti Bessel/Hankels functions.

-----------------------------------------

nonlinear coupled-mode equations solve: NL_K_Thm<-IO_K_Thm->para

IO_K_Thm: define tunable parameters, estimate cut-off power, resonance power, can sweep frequency, input power and coupling coefficient 

NL_K_Thm: function of NL coupled-mode equations.
This is to express the coupled complex nonlinear mode equations into 4 real equations, see my notebook ->NL_K_Thm

------------------------------------

mixed program to estimate quantities in Fig3 c-d: change_l->detuning

change_l: a mixture of para and IO_K_Thm, calculate for different l, only use the on resonance solution.

---------------------------------------

save_samp->detuning
to save sample given l and R.

---------------------------------------

matrices:
sample=30_166, cannot be used, delt>0
sample2=30_167
sample3=30_168
sampR_l1, as the name

------------------------------------------------------------------------------------------------
old programs:

3. kappa_cal£ºgiven l, p, n, R, calculate coupling coefficient without detuning.
Now it is not used because of different dimension used in PRA2008 and my work.

4. NL_eq: I deduced the coupled mode equations to a pair of complex nonlinear equations. Using steady state analysis, they can be further simplified into one complex nonlinear equation -> two real nonlinear equations. (see the photo 'nonlinearCoupledModeEqs')

5. IO_nonlinear: NL_eq solver (fsolve), yields Pout