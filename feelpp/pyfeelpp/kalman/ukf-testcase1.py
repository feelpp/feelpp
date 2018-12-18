def run():
    import numpy as np
    from numpy import sum, mean, sqrt
    import matplotlib.pyplot as plt
    from UKF import Filter

    def f(X,time,dt):
        Y = X + np.random.normal(0,dt)
        return Y
    
    def h(X):                                          
        Y = X #+ np.random.normal(0,dt)
        return Y

    dt = 1e-01
    Time = dt*np.arange(0,63.8)
    Signal = np.sin(Time) # + np.random.normal(0,0.1,3600) + 2
    
    F = Filter
    F.set(F,1,1,0.5,dt,f,h)
    F.readsignal(F,Signal)
    F.filter(F)

    print("L2 relative error :",np.linalg.norm(Signal-F.forecast)/np.linalg.norm(F.forecast))

    plt.plot(Time,F.forecast[0,:],label="Predicted trajectory")
    plt.plot(Time,Signal,label="Signal")
#    plt.plot(Time[0:-2],np.sin(theta*Time[0:-2])+2,label="True trajectory")
    leg = plt.legend(loc='upper right')
    leg.get_frame().set_alpha(0.5)
    plt.show()

    from matplotlib2tikz import save as tikz_save
    
#    tikz_save("ukftc1.tex")

#    plt.plot(np.log(Time[0:-2]),F.X[1,0:-2],label="R1")
#    plt.plot(np.log(1+Time[0:-2]),F.X[2,0:-2],label="R2")
#    plt.plot(np.log(1+Time[0:-2]),F.X[3,0:-2],label="C")
#    leg = plt.legend(loc='upper right')
#    leg.get_frame().set_alpha(1)
#    plt.show()

run()
