import numpy as np
from numpy import pi, sqrt




def coeff(et,eta):
   
    a0=(et**2-1)**(1/2)
    a1=-1/6*(7*et**2*eta-6*et**2-eta)/(et**2-1)**(1/2)
    a2=(1/72*(365*et**4*eta**2-1539*et**4*eta+1242*et**4-58*et**2*eta**2-690*et**2*eta
        -792*et**2+17*eta**2+69*eta+738)/(et**2-1)**(3/2))
    a3=(369/64/(et**2-1)**(5/2)*((-68204/29889*et**6-15004/9963*et**4-1580/9963*
        et**2-412/29889)*eta**3+(34700/3321*et**6+27196/3321*et**4+3412/3321*et**2+17924/
        3321)*eta**2+(-19700/1107*et**6-4036/1107*et**4+(-1772228/38745+pi**2)*et**2+1/3*pi**
        2-13916/1107)*eta+5264/1107+14096/1107*et**6-1456/123*et**4+7312/369*et**2))
    return (a0,a1,a2,a3)
    




def cal_xQ(b1,a0,a1):
    b=b1
    return (a0/(b-a1))


# In[5]:


def cal_x(x,b1,a0,a1,a2,a3,order):
    
    b=b1
    if order<=1:
        return a0/(b-a1)
    if order==2:
        return a0/(b-a1-a2*x)
    if order==3:
        return a0/(b-a1-a2*x-a3*x**2)


# In[26]:


def get_x(et,eta,b1,order):
    
    b=b1
    a0, a1, a2, a3 =coeff(et,eta)
    #a1=coeff(et,eta)[1]
    #a2=coeff(et,eta)[2]
    #a3=coeff(et,eta)[3]
    x0=cal_xQ(b1,a0,a1)
    tol = 1e-15
    diff = 1
    x1 = x0
    step = 0
    while diff>tol:
        x2 = cal_x(x1,b1,a0,a1,a2,a3,order)
        diff = x2-x1
        x1 = x2
        step += 1
    return (x1, step)
    
    
def get_x_PN(et,eta,b):
    A1=(1/6*(-7*et**2*eta+6*et**2+eta)/b/(et**2-1)**(1/2))
    A2=(1/72*(463*et**4*eta**2-1707*et**4*eta+1314*et**4-86*et**2*eta**2-666*et**2*eta-792*et**
        2+19*eta**2+69*eta+738)/b**2/(et**2-1))
    A3=(1/181440/b**3/(et**2-1)**(3/2)*(5894560*et**6*eta**3-28004760*et**6*eta**2+41840820*et**6
        *eta+480480*et**4*eta**3-22891680*et**6-12063240*et**4*eta**2-1046115*pi**2*et**2*
        eta+389340*et**4*eta+406560*et**2*eta**3+18370800*et**4+259560*et**2*eta**2-348705*pi**2
        *eta+54835596*et**2*eta-7840*eta**3-26308800*et**2-5733000*eta**2+12220740*eta-
        4974480))
    FAC=sqrt(et**2-1)/b
    return FAC*(1+A1+A2+A3)




def get_b(et,eta,x):
    a0=coeff(M1,et,eta)[0]
    a1=coeff(M1,et,eta)[1]
    a2=coeff(M1,et,eta)[2]
    a3=coeff(M1,et,eta)[3]
    B=a0/x+a1+a2*x+a3*x**2
    return(B)
    

# for i in range(0,4):
#     print(get_x(1.1,0.2,100,i))
