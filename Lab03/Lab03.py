import matplotlib.pyplot as plt
import numpy as np

def metodaTrapezow(xn,vn,deltaT,alfa):
    gamma=10**-10
    xn1=xn
    vn1=vn
    g=alfa*(1-xn**2)*vn-xn
    f=vn
    while True:
        f1=vn1
        g1=alfa*(1-xn1**2)*vn1-xn1
        F=xn1-xn-(deltaT/2)*(f+f1)
        G=vn1-vn-(deltaT/2)*(g+g1)

        a11=1
        a12=(-deltaT)/2
        a21=((-deltaT)/2)*((-2)*alfa*xn1*vn1-1)
        a22=1-(deltaT/2)*alfa*(1-(xn1**2))
        
        deltaX=((-F)*a22-(-G)*a12)/(a11*a22-a12*a21)
        deltaV=(a11*(-G)-a21*(-F))/(a11*a22-a12*a21)

        xn1=deltaX+xn1
        vn1=deltaV+vn1

        if abs(deltaV)<gamma and abs(deltaX)<gamma :
            break
    ListXV=[xn1,vn1]
    return ListXV
    
def metodaRK2(xn,vn,deltaT,alfa,znak):
    k1x=vn
    k1v=alfa*(1-xn**2)*vn-xn
    
    k2x=vn+deltaT*k1v
    k2v=alfa*(1-(xn+deltaT*k1x)**2)*(vn+deltaT*k1v)-(xn+deltaT*k1x)
    if znak=='x':
        return xn+(deltaT/2)*(k1x+k2x)
    else:
        return vn+(deltaT/2)*(k1v+k2v)

#Trapezy rozwiazanie
def trapezRozw(x0,v0,deltaT,S,p,tmax,alfa,TOL,tabT,tabDeltaT,tabX,tabV):
    xn1=0
    vn1=0
    xn2=0
    vn2=0
    t=0
    while True:
        #dwa kroki czasowe deltaT0
        XY=metodaTrapezow(x0,v0,deltaT,alfa)
        xn1=XY[0]
        vn1=XY[1]
        XY=metodaTrapezow(xn1,vn1,deltaT,alfa)
        xn1=XY[0]
        vn1=XY[1]
        #jeden krok 2*deltaT0
        XY=metodaTrapezow(x0,v0,2*deltaT,alfa)
        xn2=XY[0]
        vn2=XY[1]
        #Liczenie Ex i Ev
        Ex=(xn2-xn1)/(2**p-1)
        Ev=(vn2-vn1)/(2**p-1)
        #print(str(Ex)+" "+str(Ev))
        #Sprawdzanie akceptowalnosci wyniku
        if (max(abs(Ex),abs(Ev))<TOL):
            t=t+2*deltaT
            x0=xn2
            v0=vn2
            #print(str(t)+" "+str(deltaT)+" "+str(x0)+" "+str(v0))
            tabT.append(t)
            tabDeltaT.append(deltaT)
            tabX.append(x0)
            tabV.append(v0)
        #Zmiana kroku czasowego
        deltaT=(((S*TOL)/(max(abs(Ex),abs(Ev))))**(1/(p+1)))*deltaT
        #Warunek stopu
        if t>=tmax:
            break
#RK2 rozwiazanie
def RK2Rozw(x0,v0,deltaT,S,p,tmax,alfa,TOL,tabT,tabDeltaT,tabX,tabV):
    xn1=0
    vn1=0
    xn2=0
    vn2=0
    t=0
    while True:
        #dwa kroki deltaT
        xn1=metodaRK2(x0,v0,deltaT,alfa,'x')
        vn1=metodaRK2(x0,v0,deltaT,alfa,'v')
        xn1=metodaRK2(xn1,vn1,deltaT,alfa,'x')
        vn1=metodaRK2(xn1,vn1,deltaT,alfa,'v')
        #print(str(xn1)+" "+)
        #raz krok 2*deltaT
        xn2=metodaRK2(x0,v0,2*deltaT,alfa,'x')
        vn2=metodaRK2(x0,v0,2*deltaT,alfa,'v')
        #Liczenie Ex i Ev
        Ex=(xn2-xn1)/(2**p-1)
        Ev=(vn2-vn1)/(2**p-1)
        #print(str(Ex)+" "+str(Ev))
        if(max(abs(Ex),abs(Ev))<TOL):
            t=t+2*deltaT
            x0=xn2
            v0=vn2
            #print(str(t)+" "+str(deltaT)+" "+str(x0)+" "+str(v0))
            tabT.append(t)
            tabDeltaT.append(deltaT)
            tabX.append(x0)
            tabV.append(v0)
        #Zmiana kroku czasowego
        deltaT=(((S*TOL)/(max(abs(Ex),abs(Ev))))**(1/(p+1)))*deltaT
        #Warunek stopu
        if t>=tmax:
            break

x0=0.01
v0=0
deltaT0=1
S=0.75
p=2
tmax=40
alfa=5
TOL=10**-2
#Listy przechowuja wyniki dla TOL=10**-2
tabX2=[]
tabV2=[]
tabDeltaT2=[]
tabT2=[]
RK2Rozw(x0,v0,deltaT0,S,p,tmax,alfa,TOL,tabT2,tabDeltaT2,tabX2,tabV2)

TOL=10**-5
#Listy przechpwujace wyniki dla TOL=10**-5
tabX5=[]
tabV5=[]
tabDeltaT5=[]
tabT5=[]
RK2Rozw(x0,v0,deltaT0,S,p,tmax,alfa,TOL,tabT5,tabDeltaT5,tabX5,tabV5)
#Wykresy dla RK2
#x(t)
plt.figure(1)
#TOL=10**-2
plt.plot(tabT2,tabX2)
#TOL=10**-5
plt.plot(tabT5,tabX5)
plt.xlabel("t")
plt.ylabel("x")
plt.title("Wykres x(t)")

#v(t)
plt.figure(2)
#TOL=10**-2
plt.plot(tabT2,tabV2)
#TOL=10**-5
plt.plot(tabT5,tabV5)
plt.xlabel("v")
plt.ylabel("x")
plt.title("Wykres v(t)")

#deltaT(t)
plt.figure(3)
#TOL=10**-2
plt.plot(tabT2,tabDeltaT2)
#TOL=10**-5
plt.plot(tabT5,tabDeltaT5)
plt.xlabel("deltaT")
plt.ylabel("x")
plt.title("Wykres deltaT(t)")

#v(x)
plt.figure(4)
#TOL=10**-2
plt.plot(tabX2,tabV2)
#TOL=10**-5
plt.plot(tabX5,tabV5)
plt.xlabel("x")
plt.ylabel("v")
plt.title("Wykres v(x)")
#plt.show()
#Rozwiazanie Trapezow
#TOL=10**-2
TOL=10**-2
tabX2=[]
tabV2=[]
tabDeltaT2=[]
tabT2=[]
trapezRozw(x0,v0,deltaT0,S,p,tmax,alfa,TOL,tabT2,tabDeltaT2,tabX2,tabV2)
TOL=10**-5
tabX5=[]
tabV5=[]
tabDeltaT5=[]
tabT5=[]
trapezRozw(x0,v0,deltaT0,S,p,tmax,alfa,TOL,tabT5,tabDeltaT5,tabX5,tabV5)
#Wykresy rozwiazanie trapezow
#x(t)
plt.figure(5)
#TOL=10**-2
plt.plot(tabT2,tabX2)
#TOL=10**-5
plt.plot(tabT5,tabX5)
plt.xlabel("t")
plt.ylabel("x")
plt.title("Wykres x(t)")
#v(t)
plt.figure(6)
#TOL=10**-2
plt.plot(tabT2,tabV2)
#TOL=10**-5
plt.plot(tabT5,tabV5)
plt.xlabel("v")
plt.ylabel("x")
plt.title("Wykres v(t)")
#deltaT(t)
plt.figure(7)
#TOL=10**-2
plt.plot(tabT2,tabDeltaT2)
#TOL=10**-5
plt.plot(tabT5,tabDeltaT5)
plt.xlabel("deltaT")
plt.ylabel("x")
plt.title("Wykres deltaT(t)")
#v(x)
plt.figure(8)
#TOL=10**-2
plt.plot(tabX2,tabV2)
#TOL=10**-5
plt.plot(tabX5,tabV5)
plt.xlabel("x")
plt.ylabel("v")
plt.title("Wykres v(x)")
plt.show()
