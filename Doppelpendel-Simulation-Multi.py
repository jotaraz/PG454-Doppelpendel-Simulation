import numpy as np
import matplotlib.pyplot as plt
import imageio

#from PIL import Image, ImageDraw

Phi1 = []
Phi2 = []
Omega1 = []
Omega2 = []

dt = 1e-4 #Zeitintervall in s
g = 9.81 #Ortsfaktor in m/s^2
R1 = 0.35 #Radius des ersten Arms, in m
R2 = 0.35 #Radius des zweiten Arms, in m
m1 = 0.5 #Masse am mittleren Gelenk, in kg
m2 = 0.5 #Masse am unteren Gelenk, in kg

zp1 = 3.0 #Zeitpunkte der einzelnen Aufnahmen
zp2 = 7.9
zp3 = 11.0

xs = 400 #Groeße des GIFs in Pixeln
ys = 400

epsilon = 1e-9 #Anfangsabweichung in Phi1 in rad

Gesamt_Zeit = 5 #Die Simulation laeuft bis t = 50 s
Length = int(Gesamt_Zeit/dt) #Anzahl an Iterationen bis zum Ende
count = 0

#gif = []
f = 280.0 #Umrechnung von m in px, fuer das GIF

zeit_schweif = 0.1 #Die Zeit, die durch den Schweif abgedeckt wird

i = 0
n_tail = 50 #Anzahl an Punkten, die den Schweif bilden
dist_tail_p = int(zeit_schweif/(n_tail*dt)) #Zeitlicher Abstand zwischen zwei Punkten des Schweifs

zeit_bild = 5e-2 #Alle wie viel (virtuelle) Sekunden wird ein neues Bild gezeigt, kann auch 1 sein, ist aber dann rechenaufwendiger
zahl = int(zeit_bild/dt) #Iterationen, zwischen zwei Bildern

def calculate(p1, p2, w1, w2): #berechnet die Winkelbeschleunigungen auf Basis der Winkel p1, p2 und der Winkelgeschwindigkeiten w1, w2
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
     for i in range(nx): #ordnet jedem Pixel (entlang der x-Achse) eine y-Koordinate zu
          X.append([int(xmin + i), int(ymin + i*t)])
     if(y1 < y2):
          xmin = x1
          ymin = y1
     else:
          xmin = x2
          ymin = y2
     for i in range(ny): #ordnet jedem Pixel (entlang der y-Achse) eine x-Koordinate zu
          X.append([int(xmin + i/t), int(ymin + i)])
     return X
          

def create_gif(I1, I2): #Erstellt ein GIF, was das Doppelpendel darstellt, mit den Winkeln I1 des ersten Gelenkes und den Winkeln I2 des zweiten Gelenkes
     gif = [] #Array, was die Bilder speichert
     i = 0
     while(i < len(I1)):
          data = np.empty([xs, ys])*0 #leere Matrix, ist ein schwarzer Frame
          for k in range(N_Pendel): #geht alle Pendel (z.B. 3) durch
               A = draw_dots(I1[i][k], I2[i][k]) 
               data[xs//2, ys//2] = 255 #Mittelpunkt des Bildes (Aufhaengungspunkt des Pendels) wird weiß markiert
               data[xs//2-int(f*A[1][1]), ys//2+int(f*A[0][1])] = 255 #Gelenke werden weiß markiert
               data[xs//2-int(f*A[1][2]), ys//2+int(f*A[0][2])] = 255
               
               L = line(xs//2, xs//2-int(f*A[1][1]), ys//2, ys//2+int(f*A[0][1])) #verbindende Linie zwischen Mittelpunkt und erstem Gelenk
               for j in range(len(L)):
                    data[L[j][0], L[j][1]] = 255
               L = line(xs//2-int(f*A[1][1]), xs//2-int(f*A[1][2]), ys//2+int(f*A[0][1]), ys//2+int(f*A[0][2])) #verbindende Linie zwischen erstem und zweitem Gelenk
               for j in range(len(L)):
                    data[L[j][0], L[j][1]] = 255
               if(i > n_tail*dist_tail_p): #zeichnet Schweif
                    j = 0
                    while(j < n_tail*dist_tail_p):
                         A = draw_dots(I1[i-j][k], I2[i-j][k])
                         data[xs//2-int(f*A[1][1]), ys//2+int(f*A[0][1])] = int(230-230*j/(n_tail*dist_tail_p))
                         data[xs//2-int(f*A[1][2]), ys//2+int(f*A[0][2])] = int(230-230*j/(n_tail*dist_tail_p))
                         j += dist_tail_p
                         
          gif.append(data) #haengt den Frame an das Array an
          i += zahl
     imageio.mimsave('Doppelpendel-multi-'+str(Gesamt_Zeit)+'.gif', gif) #Erstellt ein GIF

N_Pendel = 3 #Es werden 3 Pendel simuliert

#Startbedingungen
Phi1.append(np.array([np.pi*3/4, np.pi*3/4+epsilon, np.pi*3/4+2*epsilon]))
Phi2.append(np.array([np.pi*3/4, np.pi*3/4+epsilon, np.pi*3/4+2*epsilon]))
Omega1.append(np.array([0, 0, 0]))
Omega2.append(np.array([0, 0, 0]))

#while(count < Length):
for count in range(Length): #Berechnet 'Length' Iterationen der Doppelpendel
     a = calculate(Phi1[-1], Phi2[-1], Omega1[-1], Omega2[-1])
     Omega1.append(Omega1[-1] + dt*a[0])
     Omega2.append(Omega2[-1] + dt*a[1])
     Phi1.append(Phi1[-1] + dt*Omega1[-1])
     Phi2.append(Phi2[-1] + dt*Omega2[-1])
     #count += 1
    
     if(count % zahl == 0): #Alle zahl Iterationen wird der Graph neu gezeichnet
        plt.plot([-R1-R2, R1+R2], [0, 0], 'k', label='$t$ = '+str(count*dt)+' s')
        plt.plot([0, 0], [0, -R1-R2], 'k')
        
        if(count > n_tail*dist_tail_p):
            j = 0
            while(j < n_tail*dist_tail_p):
                A = draw_dots(Phi1[count-j][0], Phi2[count-j][0])
                plt.plot(A[0], A[1], '.k')
                A = draw_dots(Phi1[count-j][1], Phi2[count-j][1])
                plt.plot(A[0], A[1], '.k')
                A = draw_dots(Phi1[count-j][2], Phi2[count-j][2])
                plt.plot(A[0], A[1], '.k')
                j += dist_tail_p
        
        A0 = draw_dots(Phi1[-1][0], Phi2[-1][0])
        A1 = draw_dots(Phi1[-1][1], Phi2[-1][1])
        A2 = draw_dots(Phi1[-1][2], Phi2[-1][2])
        
        plt.plot(A0[0], A0[1], label='Pendel Nummer 1')
        plt.plot(A0[0], A0[1], '.')
        plt.plot(A1[0], A1[1], label='Pendel Nummer 2')
        plt.plot(A1[0], A1[1], '.')
        plt.plot(A2[0], A2[1], label='Pendel Nummer 3')
        plt.plot(A2[0], A2[1], '.')

        plt.axis('equal')
        plt.legend()

        plt.draw()

        if(count*dt == zp1):
             plt.savefig('Dp-'+str(zp1)+'.pdf')
        elif(count*dt == zp2):
             plt.savefig('Dp-'+str(zp2)+'.pdf')
        elif(count*dt == zp3):
             plt.savefig('Dp-'+str(zp3)+'.pdf')
        
        plt.pause(1e-39)
        plt.cla()
        

create_gif(Phi1, Phi2)

