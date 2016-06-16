import numpy as np
from SJC_Plot import PlotLine1
FileNameX = "PtsSol.txt"
FileNameSolW = "SolW.txt"
FileNameSolU = "SolU.txt"
X = np.loadtxt(FileNameX)
SolW = np.loadtxt(FileNameSolW)
SolU = np.loadtxt(FileNameSolU)
FigName = "rho.png"
PlotLine1(X, SolW[:, 0], FigName)
FigName = "vel.png"
PlotLine1(X, SolW[:, 1], FigName)
FigName = "pressure.png"
PlotLine1(X, SolW[:, 2], FigName)

