import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
import os
import csv

def CreateDirectory(foldername):
    try: #primero se intenta guardar con un nombre sencillo
        os.mkdir(foldername)
        subfolder = foldername
    except:
        filepend = True
        trycount = 1
        while filepend:
            try: #se añade una cola al nombre de carpeta si la carpeta ya existe
                os.mkdir(foldername+' ('+str(trycount)+')')
                filepend = False
                subfolder = foldername+' ('+str(trycount)+')'
            except:
                trycount += 1
    return subfolder #el nombre de la carpeta creada

def InitializeLattice(size):
    startmode = int(entMode.get())
    #La matriz lattice contiene los espines, rodeados de filas y columnas de ceros
    if startmode==1:
        lattice = np.ones((size+2,size+2))
    elif startmode ==2:
        lattice = np.ones((size+2,size+2))
        lattice[:,:] = -1
    else:
        lattice = np.ones((size+2,size+2))
        for i in range(size):
            for j in range(size):
                chooser = np.random.random()
                lattice[i+1][j+1] -= 2*(chooser<0.5)
    #se rodea la malla con los valores de los lados opuestos (condición toroidal)
    lattice[:,0] = lattice[:,size]
    lattice[0,:] = lattice[size,:]
    lattice[size+1,:] = lattice[1,:]
    lattice[:,size+1] = lattice[:,1]
    return lattice

def PlotAndSaveMatrix(isingmatrix,subdirectory,iteration):
    size = len(isingmatrix[0])-2
    plotX = []
    plotY = []
    for x in range(size):
        for y in range(size):
            if isingmatrix[x+1,y+1]>0: #se confirma si el punto representa un +1
                plotX.append(x)
                plotY.append(y)
    endPath = subdirectory+'/'+'Ising '+str(iteration)+'.png'
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(plotX, plotY,'b*')
    ax.set_xlim([0,size-1])
    ax.set_ylim([0,size-1])
    ax.set_title('n='+str(iteration))
    fig.savefig(endPath)
    plt.close('all')
    return

def Grapher(graphMatrix,xTitle,yTitle,plotLabel,subdirectory):
    endPath = subdirectory+'/'+plotLabel+'.png'
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(graphMatrix[0], graphMatrix[1])
    ax.set_xlabel(xTitle)
    ax.set_ylabel(yTitle)
    ax.set_title(plotLabel)
    fig.savefig(endPath)
    plt.close('all')
    return

def Writer(savematrix,writepath):
    initialrow = ['n','E','M']
    rownumber = len(savematrix[0])
    colnumber = len(savematrix)
    spampath = writepath+'/results.csv'
    with open(spampath, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(initialrow)
        for row in range(rownumber):
            currentrow = []
            for col in range(colnumber):
                currentrow.append(savematrix[col][row])
            spamwriter.writerow(currentrow)
    return

def CalculateM(isingmatrix):
    size = len(isingmatrix[0])-2
    Mtotal = 0
    for i in range(size): #los 'bordes' en 0, size+1 de la matriz no aportan a M
        for j in range(size):
            Mtotal += isingmatrix[i+1][j+1] #se suman los espines totales
    return Mtotal

def CalculateE(isingmatrix,J):
    size = len(isingmatrix[0])-2
    Etotal = 0
    for m in range(size):
        for n in range(size):
            Etotal += isingmatrix[n+1][m+1]*isingmatrix[n+2][m+1]+isingmatrix[m+1][n+1]*isingmatrix[m+1][n+2]
    Etotal *= -J
    return Etotal

def CalculateResultSingleSim():
    #se obtienen los valores importantes de la interfaz
    rngseed = entrngSeed.get()
    size = int(entSize.get()) #tamaño de la malla a modificar
    cycles = int(entCycles.get()) #cantidad de iteraciones de la simulación
    frameskipPrev = entFrameskip.get() #cada cuántas iteraciones se debe graficar el modelo
    skipsize = int(entEandMskip.get()) #cada cuántas iteraciones se guardan los valores de E y M
    J = int(entJ.get())
    kbT = int(entkbT.get())
    decideGraph = bool(entDecideGraph.get())
    #se procesan los valores de la interfaz
    if rngseed != '': #se toma el valor de seed para el generador rng
        np.random.seed(int(rngseed))
    if frameskipPrev=='': #se mide si se va a graficar el modelo 2D
        frameskip = 1
        decideGraph = False
    else:
        frameskip = int(frameskipPrev)
    if decideGraph: #si se quiere graficar, se crea el directorio para hacerlo
        graphSavePath = CreateDirectory('Images')
    resultSavePath = CreateDirectory('Results') #se crea el directorio para almacenar M y E
    #se desarrolla el modelo de Ising propiamente
    lattice = InitializeLattice(size)
    E = []
    M = []
    N = []
    for n in range(cycles):
        if n%skipsize==0: #se almacenan los valores relevantes de E y M
            E.append(CalculateE(lattice,J))
            M.append(CalculateM(lattice))
            N.append(n)
        if n%frameskip==0 and decideGraph: #se decide si graficar el estado actual
            PlotAndSaveMatrix(lattice,graphSavePath,n)
        x, y, rngchoose = np.random.random(3) #se escoge el punto aleatorio a modificar
        i, j = int(x*(size-1)), int(y*(size-1)) #el punto [i+1, j+1] es el seleccionado
        deltaE = 2*J*(lattice[i,j+1]+lattice[i+2,j+1]+lattice[i+1,j]+lattice[i+1,j+2])
        Rchoose = np.exp(-deltaE/kbT) #el valor del algoritmo de Metrópolis para comparar
        lattice[i+1][j+1] *= (-1)**(rngchoose<=Rchoose)
        if (i == 0) or (i==(size-1)): #se aplica la condición toroidal de frontera en filas
            lattice[0][j+1] = lattice[size][j+1]
            lattice[size+1][j+1] = lattice[1][j+1]
        if (j == 0) or (j==(size-1)): #se aplica la condición toroidal de frontera en columnas
            lattice[i+1][0] = lattice[i+1][size]
            lattice[i+1][size+1] = lattice[i+1][1]
    #se almacenan los valores finales:
    if decideGraph:
        PlotAndSaveMatrix(lattice,graphSavePath,cycles)
    E[len(E)-1] = CalculateE(lattice,J)
    M[len(M)-1] = CalculateM(lattice)
    N[len(N)-1] = cycles
    Grapher([N,E],'n','E','E vs n',resultSavePath)
    Grapher([N,M],'n','M','M vs n',resultSavePath)
    Writer([N,E,M],resultSavePath)
    return

#EN ADELANTE SE PROGRAMA SOLAMENTE LA INTERFAZ DE USUARIO

window = tk.Tk() #crear ventana
window.title('Simulador de modelo de Ising en 2D') #título de ventana

#crear zona de entradas
frmEntry = tk.Frame(master=window)
frmEntry.grid(row=0, column=0, padx=10)

#crear entrada de size
lblSize = tk.Label(master=frmEntry, text="Tamaño de la malla:")
entSize = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblSize.grid(row=0, column=0, pady=10, padx=20, sticky="w")
entSize.grid(row=0, column=1, pady=10, padx=20, sticky="w")

#crear entrada de cycles
lblCycles = tk.Label(master=frmEntry, text="Número de ciclos:")
entCycles = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblCycles.grid(row=1, column=0, pady=10, padx=20, sticky="w")
entCycles.grid(row=1, column=1, pady=10, padx=20, sticky="w")

#crear entrada de J
lblJ = tk.Label(master=frmEntry, text="J:")
entJ = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblJ.grid(row=2, column=0, pady=10, padx=20, sticky="w")
entJ.grid(row=2, column=1, pady=10, padx=20, sticky="w")

#crear entrada de kbT
lblkbT = tk.Label(master=frmEntry, text="kbT:")
entkbT = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblkbT.grid(row=3, column=0, pady=10, padx=20, sticky="w")
entkbT.grid(row=3, column=1, pady=10, padx=20, sticky="w")

#crear entrada de radiobutton
entMode = tk.IntVar()
frmMode = tk.Frame(master=frmEntry)
lblMode = tk.Label(master=frmMode, text="Modo de inicialización de la matriz:")
R1 = tk.Radiobutton(master=frmMode, text="Matriz de 1", variable=entMode, value=1)
R2 = tk.Radiobutton(master=frmMode, text="Matriz de -1", variable=entMode, value=2)
R3 = tk.Radiobutton(master=frmMode, text="Matriz aleatoria", variable=entMode, value=3)
#organizar en el grid
frmMode.grid(row=4, column=0, padx=10, pady=10)
lblMode.grid(row=0, column=0, pady=10, padx=20, sticky="w")
R1.grid(row=0, column=1, pady=10, padx=20, sticky="w")
R2.grid(row=1, column=1, pady=10, padx=20, sticky="w")
R3.grid(row=2, column=1, pady=10, padx=20, sticky="w")

#crear entrada de EandMskip
lblEandMskip = tk.Label(master=frmEntry, text="Iteraciones entre cada valor de E y M almacenado:")
entEandMskip = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblEandMskip.grid(row=5, column=0, pady=10, padx=20, sticky="w")
entEandMskip.grid(row=5, column=1, pady=10, padx=20, sticky="w")

#crear entrada de rngseed
lblrngSeed = tk.Label(master=frmEntry, text="Seed del generador aleatorio de números (opcional):")
entrngSeed = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblrngSeed.grid(row=6, column=0, pady=10, padx=20, sticky="w")
entrngSeed.grid(row=6, column=1, pady=10, padx=20, sticky="w")

#crear entrada de graficación (checkbutton)
entDecideGraph = tk.IntVar()
cbGraph = tk.Checkbutton(master=frmEntry, text="Almacenar imágenes", variable=entDecideGraph, onvalue=1, offvalue=0, height=5, width=20)
#organizar en el grid
cbGraph.grid(row=7, column=0, pady=10, padx=20, sticky="w")

#crear entrada de frameskip
lblFrameskip = tk.Label(master=frmEntry, text="Iteraciones entre cada imagen:")
entFrameskip = tk.Entry(master=frmEntry, width=10)
#organizar en el grid
lblFrameskip.grid(row=8, column=0, pady=10, padx=20, sticky="w")
entFrameskip.grid(row=8, column=1, pady=10, padx=20, sticky="w")

#crear botón de inicialización
btnRunMain = tk.Button(master=frmEntry, text="Simular", command=CalculateResultSingleSim)
#organizar en el grid
btnRunMain.grid(row=9, column=1, sticky="e", pady=20, padx=20)

#correr la aplicación
window.mainloop()

