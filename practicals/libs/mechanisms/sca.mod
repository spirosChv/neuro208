TITLE High-Voltage Activated (HVA) Ca2+ current

COMMENT
Uses fixed eca instead of GHK equation.
Based on Reuveni, Friedman, Amitai and Gutnick (1993) J. Neurosci. 13: 4609-4621. doi: https://doi.org/10.1523/jneurosci.13-11-04609.1993

Changed from (AS Oct0899) ca.mod

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu
ENDCOMMENT

NEURON {
    SUFFIX sca
    USEION ca READ eca WRITE ica
    RANGE gcabar, ica
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gcabar = 0.001 (S/cm2)
    temp = 23 (degC) : original temp 
    Q10 = 2.3 (1) : temperature sensitivity
}

ASSIGNED {
    v (mV)
    celsius (degC)
    ica (mA/cm2)
    eca (mV)
    minf (1)
    hinf (1)
    mtau (ms)
    htau (ms)
    tadj (1)
}

STATE {
    m
    h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    tadj = Q10^((celsius - temp)/10(degC))
    ica = tadj*gcabar*pow(m, 2)*h*(v - eca)
} 

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

FUNCTION Exp(x) {
    if (x < -100) {
        Exp = 0
    } else {
        Exp = exp(x)
    }
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
    :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = 1(/mV)*y*(1 - x/y/2)
    } else {
        vtrap = 1(/mV)*x/(Exp(x/y) - 1)
    }
}

PROCEDURE rates(v (mV)) {  
    LOCAL alpha, beta

    : "m" activation gating variable
    alpha = 0.055(/ms)*vtrap(-(v + 35(mV)), 1(mV))
    beta = 0.94(/ms)*exp(-(v + 75(mV))/17(mV))

    mtau = 1/(alpha + beta)
    minf = alpha/(alpha + beta)

    : "h" inactivation 
    alpha = 0.000152(/ms)*exp(-(v + 13(mV))/50(mV))
    beta = 0.00217(/ms)/(exp(-(v + 15(mV))/28(mV)) + 1)

    htau = 1/(alpha + beta)
    hinf = alpha/(alpha + beta)
}