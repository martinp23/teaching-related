import symfit
import numpy as np
import matplotlib.pyplot as plt
from symfit import parameters, Parameter, Variable, variables, Fit, ODEModel, D

# GENERAL MODEL this general model takes a 2D weighted adjacency matrix, where the weights are rate constants.

def generalModel(rates,conc0,tvec=np.linspace(0, 200000, 100)):


    # make a list of parameters
    numEl = np.shape(rates)[0]

    model_dict = {}
    pp  = ()
    vv = ()


    # first create variables for initial species
    for ii in np.arange(0,numEl):
        var = Variable(chr(ii+65))
        vv = vv + (var,)

    t = Variable('t')

    kdict = {}

    # then create variables for products and rate constants (parameters)
    ik = 0
    for ii in np.arange(0,numEl):
        for ij in np.arange(0,numEl):
            if ii < ij:
                var = Variable(chr(ii+65)+chr(ij+65))
                vv = vv + (var,)
                par = Parameter('k'+chr(ii+65)+chr(ij+65),rates[ii,ij])
                pp = pp+(par,)
                # a dict so we can easily find rate constant indices later (this is hacky but if it works...)
                kdict[str(ii)+str(ij)] = ik
                kdict[str(ij)+str(ii)] = ik
                ik = ik+1

    # now create model
    ik=0
    for ii in np.arange(0,numEl):
        # this will be an expression for what's happening to the SM concentration. It's easiest if we just add each
        # relevant product forming reaction to this expression, then take its negative later on
        smexpr = 0
        for ij in np.arange(0,numEl):
            if ii<ij:
                model_dict[D(vv[numEl+ik],t)]= pp[kdict[str(ii)+str(ij)]]*vv[ii]*vv[ij]
                smexpr = smexpr + pp[kdict[str(ii)+str(ij)]]*vv[ii]*vv[ij]
                ik = ik+1
            elif ii>ij:
                # need to have this otherwise we miss a lot of contributions for B/C/D
                smexpr = smexpr + pp[kdict[str(ij)+str(ii)]]*vv[ij]*vv[ii]

        # we're still in the loop here, at the level of starting materials. This part creates d[A]/dt (etc for B, C, D)
        model_dict[D(vv[ii],t)] = -(smexpr)

    # set initial parameters: at time 0, all concentrations of products are zero and concentration of SMs is fixed 
    # (this could be changed to allow variable concs TODO)
    # while we're here, also set the arguments for the fit command later on to zero.
    initial = {t:0.0,}
    fitargs = {}
    for el in vv:
        if len(el.name) == 1:
            initial[el] = conc0
        else:
            initial[el] = 0
        fitargs[el.name] = None

    # define the model
    ode_model = ODEModel(model_dict, initial=initial)
    # honestly I don't know what this does but it seems to have no effect on results (based on my incomplete testing!)
    # it just needs to be there and not 'None'
    tdata = [0,1,2]

    # and then we fit the ODE model
    fit = Fit(ode_model,**fitargs,t=tdata)
    fit_result = fit.execute()


    # Generate some data from our fit model
    ans = ode_model(t=tvec, **fit_result.params)._asdict()

    # and plot it
    result = []
    legtxt = ()
    for ii in np.arange(numEl,len(vv),1):
        plt.plot(tvec, ans[vv[ii]], label=vv[ii].name)
        result.append(ans[vv[ii]][-1])
        legtxt = legtxt + (vv[ii].name,)
    plt.ylabel('Conc [M]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()

    resNorm = result/sum(result)
    xpos = np.arange(1,len(resNorm)+1,1)
    plt.bar(xpos,100*resNorm)
    plt.xticks(xpos,legtxt)
    plt.ylabel('%age at eq')
    plt.show()


    # enhancement, in percent, compared to equal concentrations everywhere
    resEnh = 100*((np.array(resNorm)) - 1/len(resNorm))/(1/len(resNorm))

    # rounding errors can give a spurious difference: set small values to zero
    resEnh[abs(resEnh) < 1e-5] = 0
    if (sum(abs(resEnh)) > 0):
        plt.bar(xpos,resEnh)
        plt.xticks(xpos,legtxt)
        plt.ylabel('%age at eq')
        plt.title('Enhancement / %')
        plt.show()
    else:
        print("No enhancement compared to equal rates in a fully-connected network")