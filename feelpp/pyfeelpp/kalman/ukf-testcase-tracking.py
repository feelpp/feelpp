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
    leg = plt.legend(loc='upper right')
    leg.get_frame().set_alpha(0.5)
    plt.show()

run()
