import math
pi=math.pi
def dlina(x,y,z):
    return (x**2+y**2+z**2)**0.5
def grinv(jr,jl,jw,l,w):
    l*=pi/180
    w*=pi/180
    return -math.cos(l)*math.sin(w)*Jw+math.cos(l)*math.cos(w)*jr-math.sin(l)*jl,-math.sin(w)*math.sin(l)*jw+math.sin(l)*math.cos(w)*jr+jl*math.cos(l),math.cos(w)*jw+jr*math.sin(w)
n,m=map(int,input("Введите диапазон сферических гармоник(n, m):").split())
lyambda,psi,r=map(float,input("Введите долготу, широту и расстояние от центра Земли до объекта:").split())
lyambda*=math.pi/180
psi*=math.pi/180
def dfact(x):
    if type(x) == int and x >= 0:
        result = lambda x: result(x-2) * x if x > 0 else 1
        return result(x)
P=[]
for i in range(n+2):
    P.append([])
    for j in range(m+2):
        P[i].append(0)
P[0][0]=1
for i in range(1,n,1):
    P[i][i]=(float(dfact(2*i-1)*(math.cos(pi*psi/180))**i))
for i in range(1,n+1,1):
    for j in range(0,i,1):
        if (i-2<0):
            P[i][j]=1/(i-j)*((2*i-1)*math.sin(psi*pi/180)*P[i-1][j])
        else:
            P[i][j]=1/(i-j)*((2*i-1)*math.sin(psi*pi/180)*P[i-1][j]+(i+j-1)*P[i-2][j])    
C = []
s = []
i = 1
k = 0
with open('C.txt') as file:
    while True:
        line = file.readline()
        line = line.rstrip()
        if not line:
            break
        s.append(float(line))
        k += 1
        if k == i:
            k = 0
            C.append(s)
            s = []
            i += 1
S = []
s = []
i = 1
k = 0
with open('S.txt') as file:
    while True:
        line = file.readline()
        line = line.rstrip()
        if not line:
            break
        s.append(float(line))
        k += 1
        if k == i:
            k = 0
            S.append(s)
            s = []
            i += 1
a=6378.1365
fM=6.67408*5.9742*10**13
fact=[]
fact.append(1)
for i in range(1,n+m,1):
    fact.append(i*fact[i-1])
Jr=Jw=Jl=0
for i in range(2,n,1):
    summ=summ1=summ2=0
    for j in range(0,min(i,m),1):
        if (j==0):  sigma=1
        else:   sigma=2
        rr=(2*i+1)**0.5*(sigma*fact[i-j]/fact[i+j])**0.5
        summ+=(C[i][j]*math.cos(j*pi/180*lyambda)+S[i][j]*math.sin(j*pi/180*lyambda))*P[i][j]*rr
        summ1+=(-C[i][j]*math.sin(j*lyambda*pi/180)+S[i][j]*math.cos(j*lyambda*pi/180))*j*rr*P[i][j]
        if (j+1<=i):
            Pstix=-j*math.tan(psi*pi/180)*P[i][j]*rr+(sigma*(i-j)*(i-j+1))**0.5+P[i][j]*rr
        summ2+=(C[i][j]*math.cos(j*lyambda*pi/180)+S[i][j]*math.sin(j*lyambda*pi/180))*Pstix
    Jr+=(i+1)*summ*(a/r)**i
    Jl+=summ1*(a/r)**i
    Jw+=summ2*(a/r)**i
Jr=-fM/r/r/10**6*(Jr+1)
Jl=fM/r/r/10**6/math.cos(psi*pi/180)*Jl
Jw=fM/r/10**6/r*Jw
x,y,z=map(float,grinv(Jr,Jl,Jw,lyambda,psi))
print('Составляющие:',grinv(Jr,Jl,Jw,lyambda,psi),'Модуль:',dlina(x,y,z),sep='\n')
    

        

