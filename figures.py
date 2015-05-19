import matplotlib.pyplot as plt





def plotXY(X,Y,nameX='t',nameY='y',figNum=1):
    """ Plot Y(X)
    """


    plt.figure(figNum)
    plt.plot(X,Y)
    plt.xlabel(nameX)
    plt.ylabel(nameY)
    plt.show()

    return




def plotXY_2(X1,Y1,X2, Y2, save=0, nameX1='t',nameY1='y',nameX2='t',nameY2='y',figNum=1):
    """ Plot Y1(X1) and Y2(X2) on two subfigures
    """
    
    plt.figure(figNum)
    plt.subplot(121)
    plt.plot(X1,Y1)
    plt.xlabel(nameX1)
    plt.ylabel(nameY1)
    plt.subplot(122)
    plt.plot(X2,Y2)
    plt.xlabel(nameX2)
    plt.ylabel(nameY2)
    plt.show()

    return
