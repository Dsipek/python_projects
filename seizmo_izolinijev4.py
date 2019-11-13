#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 17:38:44 2019

@author: sipek
"""
import itertools as it
import numpy as np


#unos povrsina unutar izoseite
sk=np.array([9512.384,3972.133,1599.26205,413.40939,66.2276])
for i in range(sk.size):
    sk[i]=np.sqrt(sk[i]/np.pi)
#intenziteti
Ik=np.array([4,5,6,7,8])


#definiranje konstanti
# =============================================================================
# alpha=0.005
# mu=1
# h=5
# I_0=6+0.5
# #definiranje parcijalnih derivacija
# b=1
# =============================================================================

cc=bb=bc=bd=cd=dd= 0.0
mu=np.log10(np.e)
b=1
#definiranje funkcja koje racunaju koeficijente/parcijalne derivacije
def bfun(x,h,alpha):
    return 1
def cfun(x,h,alpha):
    return -3*alpha*mu*(h/np.sqrt(h**2+x**2))-3*h*(1/np.sqrt(h**2+x**2)-np.sqrt(h**2+x**2)/h**2)/(np.log(10)*np.sqrt(h**2+x**2))
def dfun(x,h,alpha):
    return -3*mu*(np.sqrt(x**2+h**2)-h)
def afun(x,y,h,alpha,I_0):
    return I_0-3*np.log10(np.sqrt(x**2+h**2)/h)-3*mu*alpha*(np.sqrt(x**2+h**2)-h)-y

def get_vals(x):
    a=x[0][0]
    b=x[1][0]
    c=x[2][0]
    return (a,b,c)
def set_vals(delta,I_0):
    return (delta[0]+I_0[0],delta[1]+I_0[1],delta[2]+I_0[2])
#funkcija koja racuna 3 nepoznanice iz sustava, prima uređenu trojku pocetnih parametara
def popravka (v):
    I_0=v[0]
    h=v[1]
    alpha=v[2]
    cc=bb=bc=bd=cd=dd= 0.0
    for s in np.nditer(sk):
        bb+=bfun(s,h,alpha)*bfun(s,h,alpha)
        bc+=bfun(s,h,alpha)*cfun(s,h,alpha)
        cc+=cfun(s,h,alpha)*cfun(s,h,alpha)
        bd+=bfun(s,h,alpha)*dfun(s,h,alpha)
        cd+=cfun(s,h,alpha)*dfun(s,h,alpha)
        dd+=dfun(s,h,alpha)*dfun(s,h,alpha)
        ab=ac=ad=0.0    
    for x,y in it.izip(sk,Ik):
        ab+=afun(x,y,h,alpha,I_0)*bfun(x,h,alpha)
        ac+=afun(x,y,h,alpha,I_0)*cfun(x,h,alpha)
        ad+=afun(x,y,h,alpha,I_0)*dfun(x,h,alpha)    
    M=np.array([[bb,bc,bd],[bc,cc,cd],[bd,cd,dd]])
    L=np.array([-ab,-ac,-ad])
    L=L.reshape(3,1)
    X2=np.linalg.solve(M,L)
    return get_vals(X2)
#funkcija koja iterira kroz predhodnu funkciju sve dok delta_h nije manje od 1
def racunaj(D_in):
    H=popravka(D_in)
    G=set_vals(H,D_in)
    while(H[1]>0.01):
        print("a")
        H=popravka(G)
        G=set_vals(H,G)
    return G