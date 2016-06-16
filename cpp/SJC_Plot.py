import matplotlib.pyplot as plt
plt.style.use('ggplot')
# from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D

def PlotLine1(X, Y, FigName):
    fig = plt.figure()
    plt.plot(X, Y, 'o-')
    plt.title(FigName)
    fig.savefig(FigName, dpi=200)
    print 'Figure ', FigName, ' is saved.'
    fig.clf()

def Contour1(X, Y, U, FigName):
    from numpy import min, max, linspace
    fig = plt.figure()
    LevelN = 100
    XMin = min(X)
    XMax = max(X)
    YMin = min(Y)
    YMax = max(Y)
    plt.contourf(X, Y, U, LevelN)
    plt.colorbar()
    plt.title(FigName)
    # Set plot parameters
    axes = plt.gca()
    axes.set_aspect('equal')
    axes.set_xlabel("$X$")
    axes.set_ylabel("$Y$")
#     plt.set_xlim([XMin, XMax])
#     plt.set_xticks(linspace(XMin, XMax, 9))
#     plt.set_ylim([YMin, YMax])
#     plt.set_yticks(linspace(YMin, YMax, 3))
    fig.savefig(FigName, dpi=200)
    print 'Figure ', FigName, ' is saved.'
    fig.clf()

