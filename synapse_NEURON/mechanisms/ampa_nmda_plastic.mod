TITLE Excitatory Synapse with AMPA/NMDA components and Short-Term Plasticity

COMMENT
AMPA/NMDA currents: are based on double exponentials of a decay and a rise phase, with a scaling factor
	gampa, gnmda: opening/closing states in [0,1]
	tau1: rise, tau2: decay

w_ampa/w_nmda: the relevant weights. gbar_ampa and gbar_nmda can be omitted in a later version of the file

A portion of `inmda` is updating the `ica` (calcium current)

Authors: S. Chavlis, PhD (schavlis AT imbb.forth.gr)
ENDCOMMENT

NEURON {
	POINT_PROCESS AmpaNmdaP
	NONSPECIFIC_CURRENT isyn
	USEION ca WRITE ica
	RANGE tau1_ampa, tau2_ampa, eampa, w_ampa, gampa, gbar_ampa, iampa
	RANGE tau1_nmda, tau2_nmda, enmda, w_nmda, gnmda, gbar_nmda, inmda
	RANGE tauD, tauF, U0, r, u
	RANGE fca
}

UNITS {
	(uS) = (microsiemens)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	: AMPA parameters
	tau1_ampa = 1.0 (ms) : rise time constant AMPA
	tau2_ampa = 10.0 (ms) : decay time constant AMPA
	eampa = 0 (mV) : reversal potential AMPA
	w_ampa = 1 (1) : strength of the synapse AMPA
	gbar_ampa = 0 (uS) : maximum conductance AMPA
	
	: NMDA parameters
	tau1_nmda = 2.3 (ms) : rise time constant AMPA
	tau2_nmda = 95.0 (ms) : decay time constant AMPA
	enmda = 0 (mV) : reversal potential AMPA
	w_nmda = 1 (1) : strength of the synapse AMPA
	gbar_nmda = 0 (uS) : maximum conductance AMPA

	: Magnesium blockage of NMDA parameters
	gamma_mg = 0.0625 (/mV) : Magnesium Concentration factor
	eta_mg = 0.28011 (/mM) : Magnesium sensitivity of unblock, 1/3.57
	mg2o = 1.2 (mM) : extracellular magnesium concentration

	fca = 0.1 (1) : portion of inmda that is accumulated into ica

	: Short-term dynamics (Tsodkys, Markham model)
    tauD = 800 (ms)  : Depression time constant
    tauF = 800 (ms)  : Facilitation time constant
    U0 = .3  (1)   : initial `concentration` of vesicles
}

ASSIGNED {
	v (mV) : membrane voltage
	ica (nA) : total calcium current
	isyn (nA) : synaptic current that changes the membrane voltage

	gampa (1) : gating of AMPARs
	iampa (nA) : AMPA current

	gnmda (1) : gating of NMDARs
	inmda (nA) : NMDA current
}

STATE {
	A_ampa (1) : rise of AMPA
	B_ampa (1) : decay of AMPA
	A_nmda (1) : rise of NMDA
	B_nmda (1) : decay of NMDA
	r (1) : resources available
	u (1) : release probability
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gampa = B_ampa - A_ampa
	iampa = gbar_ampa*gampa*(v - eampa) : uS*mV=nA
	
	gnmda = B_nmda - A_nmda
	inmda = gbar_nmda*gnmda*(v - enmda)*mgblock(v, gamma_mg, eta_mg, mg2o) : uS*mV=nA

	isyn = iampa + inmda*(1-fca) : total synaptic current to update the membrane voltage
	ica =  fca*inmda : ica current from the NMDA synapse
}

DERIVATIVE states {
	A_ampa' = -A_ampa/tau1_ampa
	B_ampa' = -B_ampa/tau2_ampa

	A_nmda' = -A_nmda/tau1_nmda
	B_nmda' = -B_nmda/tau2_nmda

	r' = (1.0-r)/tauD
	u' = -u/tauF
}

INITIAL {
	: Check the time constants for AMPA and NMDA
	if (tau1_ampa/tau2_ampa > .9999) {
		tau1_ampa = .9999*tau2_ampa
	}
	if (tau1_ampa/tau2_ampa < 1e-9) {
		tau1_ampa = tau2_ampa*1e-9
	}

	if (tau1_nmda/tau2_nmda > .9999) {
		tau1_nmda = .9999*tau2_nmda
	}
	if (tau1_nmda/tau2_nmda < 1e-9) {
		tau1_nmda = tau2_nmda*1e-9
	}

	: set all values at zero in the start
	A_ampa = 0
	B_ampa = 0

	A_nmda = 0
	B_nmda = 0
	u = U0
	r = 1.0
}

FUNCTION factor(tau1 (ms), tau2 (ms)) {
	LOCAL tpeak, peak
	tpeak = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1) : time of peak
	peak = exp(-tpeak/tau2) - exp(-tpeak/tau1) : peak conductance
	factor = 1/peak
}

FUNCTION mgblock(v (mV), gamma (/mV), eta (/mM), mgo (mM)) {
	mgblock = 1/(1 + eta*mgo*exp(-gamma*v))
}

NET_RECEIVE(weight) {
	LOCAL r_minus, x_val

	INITIAL {
		weight = 1
	}

	:Calculate the u and R variables at each presynaptic spike

	u = u + U0*(1.0 - u)
	r_minus = r
	r = r + (-u*r)
	
	: At every pre-synaptic spike update A (rise) and B (decay) phases
	: of AMPA and NMDA
	x_val = w_ampa*weight*factor(tau1_ampa, tau2_ampa)*r_minus*u
	A_ampa = A_ampa + x_val
	B_ampa = B_ampa + x_val
    
    x_val = w_nmda*weight*factor(tau1_nmda, tau2_nmda)*r_minus*u
	A_nmda = A_nmda + x_val
	B_nmda = B_nmda + x_val
}