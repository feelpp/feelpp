def run():
    import numpy as np
    from numpy import sum, mean, sqrt
    import matplotlib.pyplot as plt
    from UKF import Filter

    def f(X,time,dt):
        Y = X
        Y[0] = theta[0]*( theta[1]*Signal[time] - theta[2]*Signal[time-1] + C*X[0]/dt )
        Y[1] += np.random.normal(0,dt)
        return Y
    
    def h(X):
        Y = X[0]/R2 + C*X[1]
        return Y

    dt = 1e-01
    Time = np.arange(0,6.4,dt)
    Signal = np.sin(Time) + np.random.normal(0,0.1,np.size(Time))
    
    R1 = 1.6e03
    R2 = 1.3e04
    C = 2.5e-05
    theta = np.array([[1/( C/dt + 1/R2 )],[1 + R1/R2 + R1*C/dt],[R1*C/dt]])

    F = Filter
    F.set(F,2,1,0.5,dt,f,h)
    F.readsignal(F,Signal)
    F.filter(F)

    realstate = 0*Time
    for i in range(1,np.size(realstate)):
        realstate[i] = theta[0]*( theta[1]*np.sin(Time[i]) - theta[2]*np.sin(Time[(i-1)]) + C*realstate[i-1]/dt )

    print("L2 relative error :",np.linalg.norm(realstate[1:]-F.X[0,:])/np.linalg.norm(F.X[0,:]))

    plt.plot(Time[1:],F.X[0,:],label="Analyzed trajectory")
    plt.plot(Time,realstate,label="Real state")
#    plt.plot(Time[0:-2],np.sin(theta*Time[0:-2])+2,label="True trajectory")
    leg = plt.legend(loc='upper right')
    leg.get_frame().set_alpha(0.5)
    plt.show()

    from matplotlib2tikz import save as tikz_save
    
    tikz_save("ukftc3.tex")

#    plt.plot(np.log(Time[0:-2]),F.X[1,0:-2],label="R1")
#    plt.plot(np.log(1+Time[0:-2]),F.X[2,0:-2],label="R2")
#    plt.plot(np.log(1+Time[0:-2]),F.X[3,0:-2],label="C")
#    leg = plt.legend(loc='upper right')
#    leg.get_frame().set_alpha(1)
#    plt.show()

run()
