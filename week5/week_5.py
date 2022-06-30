from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import sys
import numpy as np
import scipy.linalg as s

#c_ fn to combine two different coloumn and form a matrix 
def c_(x,y):
    M = np.zeros((len(x),2))
    for i in range(len(y)):
        M[i,0] = x[i]
        M[i,1] = y[i] 
    return M

#finding A nd B in y = A*exp(Bx)
def finding_B_A (error, Niter, no_of_points=0):
    log_err = np.log(error)[-no_of_points:]
    M =  c_( (np.array(range(Niter)) +1)[-no_of_points:], np.ones(log_err.shape[0]))
    x = s.lstsq(M, log_err)

    return x[0][0], x[0][1]

#Checking if external arguments provided
try:
    len(sys.argv) == 5
    Nx=int(sys.argv[1])
    Ny=int(sys.argv[2])
    radius=int(sys.argv[3])
    Niter=int(sys.argv[4])
     
except:
    print("No external arguments. Using Default parameters, Nx = 25 Ny = 25 radius = 8 Niter = 1500")
    Nx = 25; 
    Ny = 25; 
    radius = 8;
    Niter = 1500; 

phi = np.zeros((Ny,Nx))                                        #Initialising phi to 0 matrix

x = linspace(-0.5, 0.5, Nx)                                    #Definging X coords                                                                   
y = linspace(-0.5, 0.5, Ny)                                    #Definging Y coords                            
Y,X=meshgrid(y,x)                                              #Forming meshgrid

ii = where( (X *X +Y *Y )<= 0.35 * 0.35)                       #Finding 1 potential points
phi[ii]=1.0                                    

#Plotting counter of potential
figure(1)
contourf(X,Y, phi)
colorbar()
title("contour plot of the potential")
xlabel("X coords")
ylabel("Y coords")
savefig('fig1.png')
show()

#Forming eerror matrix and phi
errors = np.zeros(Niter)                                      
for k in range(Niter):
    oldphi = phi.copy()
    left_matrix = oldphi[1:-1,0:-2]
    right_matrix = oldphi[1:-1 ,2: ]
    top_matrix = oldphi[0:-2 ,1:-1 ]
    down_matrix= oldphi[2: ,1:-1 ]

    phi[1:-1, 1:-1] =0.25*(left_matrix + right_matrix + top_matrix + down_matrix )
    
    phi[:,0]=phi[:,1]
    phi[:, -1] = phi[:, -2]
    phi[ 0 ,:] = phi[1, :]
    phi[-1, :] = 0

    phi[ii] =1.0
    
    errors[k]=np.max(np.abs(phi-oldphi))

#print(phi)
Nither_as_x = np.array(range(Niter)) +1                                       #used for x coords

#error on semilog
figure(2)
title("Error on semilog")
xlabel("x coords")
ylabel("error")
grid() 
semilogy(Nither_as_x, errors)
savefig('fig2.png')
show()

#error on loglog
figure(3)
title("Error on loglog")
xlabel("x coords")
ylabel("error")
grid() 
loglog(Nither_as_x, errors, label ="real")
loglog(Nither_as_x[::50], errors[::50], "ro", label = "every 50th")
legend()
savefig('fig3.png')
show()

#fit1 and fit2
#Finding A and B for all and 500 iterations respectivly
B,A = finding_B_A(errors,Niter)
B2,A2 = finding_B_A(errors, Niter, 500)

figure(4)
title("Error on loglog with fit1 and fit2")
xlabel("x coords")
ylabel("error")
grid() 
loglog(Nither_as_x, errors, label = "true")
loglog(Nither_as_x[::100], np.exp((B*Nither_as_x + A)[::100]), "ro" , label = "fit1"  )
loglog(Nither_as_x[::100], np.exp((B2*Nither_as_x + A2)[::100]), "go", label = "fit2"   )
legend()
savefig('fig4.png')
show()

figure(5)
title("Error on semilog with fit1 and fit2 ")
xlabel("x coords")
ylabel("error")
grid() 
semilogy(Nither_as_x, errors)
semilogy(Nither_as_x[::100], np.exp((B*Nither_as_x + A)[::100]), "ro"  , label="fit1" )
semilogy(Nither_as_x[::100], np.exp((B2*Nither_as_x + A2)[::100]), "ro", label="fit2"   )
legend()
savefig('fig5.png')
show()

# Cummulative error in loglog
figure(6)
title("Cummulative max error on loglog")
xlabel("No. of iterations")
ylabel("Cummulative error")
iter = np.arange(100, 1501, 100)
loglog(iter, np.abs( (A2/B2)*np.exp(B2*(iter +0.5)) ) , "ro")
savefig('fig6.png')
show()

fig4=figure(7) # open a new figure
#title("3-D surface plot of the potential")
ax=p3.Axes3D(fig4) # Axes3D is the means to do a surface plot
ax.set_title("The 3-D surface plot of the potential")
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
xlabel("y coords")
ylabel("x coords")
ax.set_zlabel("phi")
savefig('fig7.png')
show()

# countour plot of potentials
figure(8)
title("2D Contour plot of potential")
xlabel("X")
ylabel("Y")
grid()
plot((ii[0]-Nx/2)/Nx,(ii[1]-Ny/2)/Ny,'ro')
contourf(Y,X[::-1],phi)
colorbar()
savefig('fig8.png')
show()


Jx,Jy = (1/2*(phi[1:-1,0:-2]-phi[1:-1,2:]),1/2*(phi[:-2,1:-1]-phi[2:,1:-1]))   # finding Jx and Jy

#plotting current density
figure(9)
title("Vector plot of current flow")
quiver(Y[1:-1,1:-1],-X[1:-1,1:-1],-Jx[:,::-1],-Jy)
plot((ii[0]-Nx/2)/Nx,(ii[1]-Ny/2)/Ny,'ro')
grid()
savefig('fig9.png')
show()

