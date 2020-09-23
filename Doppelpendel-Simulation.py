import numpy as np
import matplotlib.pyplot as plt
import imageio

Phi1 = []
Phi2 = []
Omega1 = []
Omega2 = []

dt = 1e-4 #Zeitintervall in s
g = 9.81 #Ortsfaktor in m/s^2
R1 = 0.351 #Radius des ersten Arms, in m
R2 = 0.24 #Radius des zweiten Arms, in m
m1 = 0.127 #Masse am mittleren Gelenk, in kg
m2 = 0.157 #Masse am unteren Gelenk, in kg

xs = 400 #Groeße des GIFs in Pixeln
ys = 400

Gesamt_Zeit = 10 #Die Simulation laeuft bis t = 10 s
Length = int(Gesamt_Zeit/dt) #Anzahl an Iterationen bis zum Ende
count = 0

f = 280.0 #Umrechnung von m in px, fuer das GIF

zeit_schweif = 0.1 #Die Zeit, die durch den Schweif abgedeckt wird

i = 0
n_tail = 50 #Anzahl an Punkten, die den Schweif bilden
dist_tail_p = int(zeit_schweif/(n_tail*dt)) #Zeitlicher Abstand zwischen zwei Punkten des Schweifs

zeit_bild = 5e-2 #Alle wie viel (virtuelle) Sekunden wird ein neues Bild gezeigt, kann auch 1 sein, ist aber dann rechenaufwendiger
zahl = int(zeit_bild/dt) #Iterationen, zwischen zwei Bildern

def calculate_new(p1, p2, w1, w2): #berechnet die Winkelbeschleunigungen auf Basis der Winkel p1, p2 und der Winkelgeschwindigkeiten w1, w2
     c1 = -m2*w1*R1*w2*R2*np.sin(p1-p2) - (m1+m2)*g*R1*np.sin(p1)
     c2 =  m2*w1*R1*w2*R2*np.sin(p1-p2) - m2*g*R2*np.sin(p2)
     a1 = (m1+m2)*R1**2
     a2 = m2*R1*R2*np.cos(p1-p2)
     a3 = -(w1-w2)*m2*R1*w2*R2*np.sin(p1-p2)
     b1 = m2*R2**2
     b2 = m2*R1*R2*np.cos(p1-p2)
     b3 = -(w1-w2)*m2*w1*R1*R2*np.sin(p1-p2)
     acc2 = (c2 - b3 - b2*c1/a1 + b2*a3/a1)/(b1 - b2*a2/a1)
     acc1 = (c1-a2*acc2-a3)/a1
     
     return [acc1, acc2]

def draw_dots(p1, p2): #berechnet die x und y Position der Gelenke aus den Winkeln des Pendels  
     x1 = R1*np.sin(p1)
     y1 = -R1*np.cos(p1)
     x2 = x1 + R2*np.sin(p2)
     y2 = y1 - R2*np.cos(p2)
     XA = [0, x1, x2]
     YA = [0, y1, y2]
     return [XA, YA]

def draw1(p1): #erstellt ein Array mit x-Koordinaten und eins mit y-Koordinaten des ersten Gelenks aus einem Array mit Winkeln
     x1 =  R1*np.sin(p1)
     y1 = -R1*np.cos(p1)
     return [x1, y1]

def draw2(p1, p2): #erstellt ein Array mit x-Koordinaten und eins mit y-Koordinaten des zweiten Gelenks aus einem Array mit Winkeln
     x1 = R1*np.sin(p1)
     y1 = -R1*np.cos(p1)
     x2 = x1 + R2*np.sin(p2)
     y2 = y1 - R2*np.cos(p2)
     return [x2, y2]

def line(x1, x2, y1, y2): #erstellt eine Linie (aus Pixeln) im GIF
     if(x1 < x2):
          xmin = x1
          ymin = y1
     else:
          xmin = x2
          ymin = y2
     nx = abs(int(x2-x1))
     ny = abs(int(y2-y1))
     X = []
     t = (y2-y1)/(x2-x1+0.0001)
     for i in range(nx):
          X.append([int(xmin + i), int(ymin + i*t)])
     if(y1 < y2):
          xmin = x1
          ymin = y1
     else:
          xmin = x2
          ymin = y2
     for i in range(ny):
          X.append([int(xmin + i/t), int(ymin + i)])
     return X
          

def create_gif(I1, I2): #Erstellt ein GIF, was das Doppelpendel darstellt, mit den Winkeln I1 des ersten Gelenkes und den Winkeln I2 des zweiten Gelenkes
     gif = []
     i = 0
     while(i < len(I1)):
     #for i in range(len(I1)):
          data = np.empty([xs, ys])*0
          A = draw_dots(I1[i], I2[i])
          data[xs//2, ys//2] = 255
          data[xs//2-int(f*A[1][1]), ys//2+int(f*A[0][1])] = 255
          data[xs//2-int(f*A[1][2]), ys//2+int(f*A[0][2])] = 255
               
          L = line(xs//2, xs//2-int(f*A[1][1]), ys//2, ys//2+int(f*A[0][1]))
          for j in range(len(L)):
               data[L[j][0], L[j][1]] = 255
          L = line(xs//2-int(f*A[1][1]), xs//2-int(f*A[1][2]), ys//2+int(f*A[0][1]), ys//2+int(f*A[0][2]))
          for j in range(len(L)):
               data[L[j][0], L[j][1]] = 255
          if(i > n_tail*dist_tail_p):
               j = 0
               while(j < n_tail*dist_tail_p):
                    A = draw_dots(I1[i-j], I2[i-j])
                    data[xs//2-int(f*A[1][1]), ys//2+int(f*A[0][1])] = int(230-230*j/(n_tail*dist_tail_p))
                    data[xs//2-int(f*A[1][2]), ys//2+int(f*A[0][2])] = int(230-230*j/(n_tail*dist_tail_p))
                    j += dist_tail_p
                    
          gif.append(data)
          i += zahl
     imageio.mimsave('Doppelpendel-Sim-'+str(Gesamt_Zeit)+'.gif', gif)

#Startbedingungen
Phi1.append(np.pi/2)
Phi2.append(0)
Omega1.append(0)
Omega2.append(0)


while(count < Length):
     a = calculate_new(Phi1[-1], Phi2[-1], Omega1[-1], Omega2[-1]) #a enthaelt die beiden Winkelbeschleunigungen
     Omega1.append(Omega1[-1] + dt*a[0]) #Berechnet die Winkelgeschwindigkeiten
     Omega2.append(Omega2[-1] + dt*a[1])
     Phi1.append(Phi1[-1] + dt*Omega1[-1]) #Berechnet die Winkel
     Phi2.append(Phi2[-1] + dt*Omega2[-1])
     count += 1

     if(count % zahl == 0): #Alle zahl iterationen wird der Graph 'geupdatet'
          plt.plot([-R1-R2, R1+R2], [0, 0], 'k', label='$t$ = '+str(count*dt)+' s')
          plt.plot([0, 0], [0, -R1-R2], 'k')
          if(count > n_tail*dist_tail_p):
               j = 0
               while(j < n_tail*dist_tail_p):
                    A = draw_dots(Phi1[count-j], Phi2[count-j])
                    plt.plot(A[0], A[1], '.k')
                    j += dist_tail_p
          
          A = draw_dots(Phi1[-1], Phi2[-1])
          plt.plot(A[0], A[1])
          plt.plot(A[0], A[1], '.')
          
          plt.axis('equal')
          plt.xlim([-R1-R2, R1+R2])
          plt.ylim([-R1-R2, R1+R2])
          plt.legend()
          plt.draw()
          plt.pause(1e-39)
          plt.cla()

create_gif(Phi1, Phi2)

print('done')
