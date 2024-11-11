import symfit
import numpy as np
import matplotlib.pyplot as plt
from symfit import parameters, Parameter, Variable, variables, Fit, ODEModel, D


# In this block, four models are individually defined. This is bad code: lots is repeated! But each is perhaps a clearer 
# of the use of symfit

# MODEL ONE this model comprises A+B--->AB (ie a simple chemical reaction)
def oneProductModel(kABval=1e-2,conc0=50e-3,tvec=np.linspace(0, 200000, 100)):
    # Here we describe a model with A+B->AB
    A, B, AB, t = variables('A, B, AB, t')
    tdata = [0,1,2]

    kAB = Parameter('kAB',kABval)  # Rate constant for formation of AB

    # here's a list of rate expressions for each component in the mixture
    model_dict = {
        D(AB,t): kAB * A * B,
        D(A,t): -kAB*A*B,
        D(B,t): -(kAB*A*B),
    }
    # here we define the ODE model and specify the start concentrations of each reagent

    ode_model = ODEModel(model_dict, initial={t: 0.0, A:conc0, B:conc0, AB:0, })

    # and then we fit the ODE model
    fit = Fit(ode_model, t=tdata, A=None, B=None, AB=None)
    fit_result = fit.execute()


    # Generate some data
    ans = ode_model(t=tvec, **fit_result.params)._asdict()
    
    # and plot it
    plt.plot(tvec, ans[AB], label='[AB]')


    #plt.plot(tvec, BCres, label='[BC]')
    #plt.scatter(tdata, adata)
    plt.ylabel('Conc [M]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()

# MODEL TWO this model comprises A+B--->AB  and B+C---->BC
def twoProductModel(kABval=1e-2,kBCval=1e-2,conc0=50e-3,tvec=np.linspace(0, 200000, 100)):
    # conc0 is initial concentration
    tdata = [0,1,2]
    # Here we describe a model with A+B->AB and B+C->BC
    A, B, C, AB, BC, t = variables('A, B, C, AB, BC, t')
    kAB = Parameter('kAB',kABval)  # Rate constant for formation of AB
    kBC = Parameter('kBC',kBCval)    # rate constant for formation of BC
    
    # here's a list of rate expressions for each component in the mixture
    model_dict = {
        D(AB,t): kAB * A * B,
        D(BC,t): kBC * B * C,
        D(A,t): -kAB*A*B,
        D(B,t): -(kAB*A*B + kBC*B*C),
        D(C,t): -kBC*B*C,

    }
    
    # here we define the ODE model and specify the start concentrations of each reagent

    ode_model = ODEModel(model_dict, initial={t: tdata[0], A:conc0, B:conc0, C:conc0, AB:0, BC:0})

    # and then we fit the ODE model
    # fit = Fit(ode_model, t=tdata, A=None, B=None, AB=None, BC=None, C=None)
    # fit_result = fit.execute()
    
    # Generate some data
    ans = ode_model(t=tvec,kAB=kABval, kBC=kBCval)._asdict()

    
    # and plot it
    plt.plot(tvec, ans[AB], label='[AB]')
    plt.plot(tvec, ans[BC], label='[BC]')
    #plt.scatter(tdata, adata)
    plt.ylabel('Conc [M]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()
    
    res = [ans[AB][-1],ans[BC][-1]]
    resNorm = res/sum(res)
    plt.bar([1,2],100*resNorm)
    plt.xticks([1,2],('[AB]','[BC]'))
    plt.ylabel('%age at eq')
    plt.show()

    # enhancement, in percent, compared to equal concentrations everywhere
    resEnh = 100*((np.array(resNorm)) - 1/len(resNorm))/(1/len(resNorm))
    # rounding errors can give a spurious difference: set small values to zero
    resEnh[abs(resEnh) < 1e-5] = 0
    if (sum(abs(resEnh)) > 0):
        plt.bar([1,2],resEnh)
        plt.xticks([1,2],('[AB]','[BC]'))
        plt.ylabel('%age at eq')
        plt.title('Enhancement / %')
        plt.show()
    else:
        print("No enhancement compared to equal rates")
    
    



# MODEL 3: this is a box where each vertex is connected to two others.
def box(kABval=1e-2,kACval=1e-2,kBDval=1e-2,kCDval=1e-2,conc0=50e-3,tvec=np.linspace(0, 200000, 100)):
    # Here we describe a model with A+B->AB
    A, B, C, Di, AB, AC, CD, BD,  t = variables('A, B, C, Di, AB, AC, CD, BD,  t')
    tdata = [0,1,2,100,1000,10000]

    kAB = Parameter('kAB',kABval)  # Rate constant for formation of AB
    kAC = Parameter('kAC',kACval)  # Rate constant for formation of AC
    kBD = Parameter('kBD',kBDval)  # Rate constant for formation of BD
    kCD = Parameter('kCD',kCDval)  # Rate constant for formation of CD



    # here's a list of rate expressions for each component in the mixture
    # here I'm calling the concentration of D as'Di' to avoid confusion
    model_dict = {
        D(AB,t): kAB * A * B,
        D(AC,t): kAC * A * C,
        D(BD,t): kBD * B * Di,
        D(CD,t): kCD * C * Di,
        D(A,t): -(kAB*A*B + kAC*A*C),
        D(B,t): -(kAB*A*B + kBD*B*Di),
        D(C,t): -(kAC*A*C+kCD*C*Di),
        D(Di,t): -(kBD*B*Di+kCD*C*Di),
    }
    # here we define the ODE model and specify the start concentrations of each reagent

    ode_model = ODEModel(model_dict, initial={t: 0.0, A:conc0, B:conc0, C:conc0, Di:conc0, AB:0,  AC:0, BD:0, CD:0,  })

    # and then we fit the ODE model
    # fit = Fit(ode_model, t=tdata, A=None, B=None, C=None, Di=None, AB=None,  AC=None, BD=None, CD=None, )
    # fit_result = fit.execute()


    # Generate some data
    ans = ode_model(t=tvec, kAB=kABval, kAC=kACval, kBD=kBDval, kCD=kCDval)._asdict()

    
    # and plot it
    plt.plot(tvec, ans[AB], label='[AB]')
    plt.plot(tvec, ans[AC], label='[AC]')
    plt.plot(tvec, ans[CD], label='[CD]')
    plt.plot(tvec, ans[BD], label='[BD]')

    plt.xlabel('Time [s]')
    plt.ylabel('Conc [M]')
    #plt.plot(tvec, BCres, label='[BC]')
    #plt.scatter(tdata, adata)
    plt.legend()
    plt.show()

    res = [ans[AB][-1],ans[AC][-1],ans[CD][-1],ans[BD][-1]]
    resNorm = res/sum(res)
    plt.bar([1,2,3,4],100*resNorm)
    plt.xticks([1,2,3,4],('[AB]','[AC]','[CD]','[BD]'))
    plt.ylabel('%age at eq')
    plt.show()

    # enhancement, in percent, compared to equal concentrations everywhere
    resEnh = 100*((np.array(resNorm)) - 1/len(resNorm))/(1/len(resNorm))
    
    # rounding errors can give a spurious difference: set small values to zero
    resEnh[abs(resEnh) < 1e-5] = 0
    if (sum(abs(resEnh)) > 0):
        yval = [1,2,3,4]
        plt.bar(yval,resEnh)
        plt.xticks(yval,('[AB]','[AC]','[CD]','[BD]'))
        plt.ylabel('%age at eq')
        plt.title('Enhancement / %')
        plt.show()
    else:
        print("No enhancement compared to equal rates")


# MODEL 4: this is a square where each vertex is connected to three others (i.e. it's fully connected)
def square(kABval=1e-2,kACval=1e-2,kBDval=1e-2,kCDval=1e-2,kBCval=1e-2,kADval=1e-2,conc0=50e-3,tvec=np.linspace(0, 200000, 100)):
    # Here we describe a model with A+B->AB
    A, B, C, Di, AB, AC, CD, BD, AD, BC, t = variables('A, B, C, Di, AB, AC, CD, BD, AD, BC, t')
    tdata = [0,1,2,100,1000,10000]

    kAB = Parameter('kAB',kABval)  # Rate constant for formation of AB
    kAC = Parameter('kAC',kACval)  # Rate constant for formation of AC
    kBD = Parameter('kBD',kBDval)  # Rate constant for formation of BD
    kCD = Parameter('kCD',kCDval)  # Rate constant for formation of CD
    kBC = Parameter('kBC',kBCval)  # Rate constant for formation of BC  ## cross-connection
    kAD = Parameter('kAD',kADval)  # Rate constant for formation of AD  ## cross-connection


    # here's a list of rate expressions for each component in the mixture
    # here I'm calling the concentration of D as'Di' to avoid confusion
    model_dict = {
        D(AB,t): kAB * A * B,
        D(AC,t): kAC * A * C,
        D(BD,t): kBD * B * Di,
        D(CD,t): kCD * C * Di,
        D(BC,t): kBC * B * C,
        D(AD,t): kAD * A * Di,
        D(A,t): -(kAB*A*B + kAC*A*C + kAD*A*Di),
        D(B,t): -(kAB*A*B + kBD*B*Di + kBC*B*C),
        D(C,t): -(kAC*A*C+kCD*C*Di+kBC*B*C),
        D(Di,t): -(kBD*B*Di+kCD*C*Di+kAD*A*Di),
    }
    # here we define the ODE model and specify the start concentrations of each reagent

    ode_model = ODEModel(model_dict, initial={t: 0.0, A:conc0, B:conc0, C:conc0, Di:conc0, AB:0, BC:0, AC:0, BD:0, CD:0, AD:0 })

    # and then we fit the ODE model
    # fit = Fit(ode_model, t=tdata, A=None, B=None, C=None, Di=None, AB=None, BC=None, AC=None, BD=None, CD=None, AD=None)
    # fit_result = fit.execute()


    # Generate some data
    ans = ode_model(t=tvec, kAB=kABval, kAC=kACval, kBD=kBDval, kCD=kCDval, kBC=kBCval, kAD=kADval)._asdict()

    
    # and plot it
    plt.plot(tvec, ans[AB], label='[AB]')
    plt.plot(tvec, ans[AC], label='[AC]')
    plt.plot(tvec, ans[CD], label='[CD]')
    plt.plot(tvec, ans[BC], label='[BC]')
    plt.plot(tvec, ans[BD], label='[BD]')
    plt.plot(tvec, ans[AD], label='[AD]')

    #plt.plot(tvec, BCres, label='[BC]')
    #plt.scatter(tdata, adata)
    plt.ylabel('Conc [M]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()

    res = [ans[AB][-1],ans[AC][-1],ans[CD][-1],ans[BC][-1],ans[BD][-1],ans[AD][-1]]
    resNorm = res/sum(res)
    plt.bar([1,2,3,4,5,6],100*resNorm)
    plt.xticks([1,2,3,4,5,6],('[AB]','[AC]','[CD]','[BC]','[BD]','[AD]'))
    plt.ylabel('%age at eq')
    plt.show()

    # enhancement, in percent, compared to equal concentrations everywhere
    resEnh = 100*((np.array(resNorm)) - 1/len(resNorm))/(1/len(resNorm))
    
    # rounding errors can give a spurious difference: set small values to zero
    resEnh[abs(resEnh) < 1e-5] = 0
    if (sum(abs(resEnh)) > 0):
        yval = [1,2,3,4,5,6]
        plt.bar(yval,resEnh)
        plt.xticks(yval,('[AB]','[AC]','[CD]','[BC]','[BD]','[AD]'))
        plt.ylabel('%age at eq')
        plt.title('Enhancement / %')
        plt.show()
    else:
        print("No enhancement compared to equal rates")