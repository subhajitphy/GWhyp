import numpy as np
from matplotlib import pyplot as plt

from numpy import sinh, cosh, tanh,  arctan, pi



from getx import get_x

# e1=1.5
# eta=0.25
# B1=700
# l0=5


def get_u_hat(l,e):

 
     alpha = (e-1)/(4*e + 0.5)

     alpha3 = alpha*alpha*alpha

     beta   = (l/2)/(4*e + 0.5)

     beta2  = beta*beta;

     z = np.zeros_like(l)

     

     z= np.cbrt(beta + np.sqrt(alpha3 + beta2))
     
     
        
     
     s=(z - alpha/z)
     s5= s*s*s*s*s
     
     ds= 0.071*s5/((1+0.45*s*s)*(1+4*s*s)*e)
     w= s+ds
        
     u = 3*np.log(w+np.sqrt(1+w*w))
     
     esu= e*np.sinh(u)
     ecu= e*np.cosh(u)
     
     fu  = -u + esu - l
     f1u = -1 + ecu  
     f2u = esu
     f3u = ecu
     f4u = esu
     f5u = ecu

     u1 = -fu/ f1u
     u2 = -fu/(f1u + f2u*u1/2)
     u3 = -fu/(f1u + f2u*u2/2 + f3u*(u2*u2)/6.0)
     u4 = -fu/(f1u + f2u*u3/2 + f3u*(u3*u3)/6.0 + f4u*(u3*u3*u3)/24.0)
     u5 = -fu/(f1u + f2u*u4/2 + f3u*(u4*u4)/6.0 + f4u*(u4*u4*u4)/24.0 + f5u*(u4*u4*u4*u4)/120.0)
     uM = (u + u5)
 
     return uM



def get_u(l,et,eta,b1,order):
    
    U=get_u_hat(l,et)
    x=get_x(et,eta,b1,order)[0]

    a0=U

    a2=(1/8*x**2*(24*(-5+2*eta)*arctan(((1+et)/(-1+et))**(1/2)*tanh(1/2*U))*(-1+et*cosh(U
        ))/(et**2-1)**(1/2)+et*(-15+eta)*eta*sinh(U))/(-1+et*cosh(U))**2)


    a3=(x**3*(1/8*et*(-4+eta)*(-60+3*(5*et**2+8)*eta-et**2*eta**2+et*(eta**2-39*eta+60)*cosh
        (U))*sinh(U)/(et**2-1)/(-1+et*cosh(U))**3-1/6720/(et**2-1)**(3/2)/(-1+et*cosh(U))*(
        et*(et**2-1)**(1/2)*(67200-3*(1435*pi**2+105*et**2-47956)*eta-105*(135*et**2+592)*
        eta**2+35*(65*et**2-8)*eta**3)*sinh(U)/(-1+et*cosh(U))+70*(8640+(123*pi**2-13184)*
        eta+960*eta**2+96*et**2*(11*eta**2-29*eta+30))*arctan(((1+et)/(-1+et))**(1/2)*tanh(
        1/2*U))+840*et**2*(et**2-1)**(1/2)*eta*(3*eta**2-49*eta+116)*(et-cosh(U))*sinh(U)/(
        -1+et*cosh(U))**2-35/2*et**3*(et**2-1)**(1/2)*eta*(13*eta**2-73*eta+23)*(-7*et**2-2+
        12*et*cosh(U)+(et**2-4)*cosh(2*U))*sinh(U)/(-1+et*cosh(U))**3)))

    if order<=1:
        return a0
    if order==2:
        return a0+a2
    if order==3:
        return a0+a2+a3


def get_u_v2(l,et,eta,x,order):
    
    U=get_u_hat(l,et)
    
    a0=U

    a2=(1/8*x**2*(24*(-5+2*eta)*arctan(((1+et)/(-1+et))**(1/2)*tanh(1/2*U))*(-1+et*cosh(U
        ))/(et**2-1)**(1/2)+et*(-15+eta)*eta*sinh(U))/(-1+et*cosh(U))**2)


    a3=(x**3*(1/8*et*(-4+eta)*(-60+3*(5*et**2+8)*eta-et**2*eta**2+et*(eta**2-39*eta+60)*cosh
        (U))*sinh(U)/(et**2-1)/(-1+et*cosh(U))**3-1/6720/(et**2-1)**(3/2)/(-1+et*cosh(U))*(
        et*(et**2-1)**(1/2)*(67200-3*(1435*pi**2+105*et**2-47956)*eta-105*(135*et**2+592)*
        eta**2+35*(65*et**2-8)*eta**3)*sinh(U)/(-1+et*cosh(U))+70*(8640+(123*pi**2-13184)*
        eta+960*eta**2+96*et**2*(11*eta**2-29*eta+30))*arctan(((1+et)/(-1+et))**(1/2)*tanh(
        1/2*U))+840*et**2*(et**2-1)**(1/2)*eta*(3*eta**2-49*eta+116)*(et-cosh(U))*sinh(U)/(
        -1+et*cosh(U))**2-35/2*et**3*(et**2-1)**(1/2)*eta*(13*eta**2-73*eta+23)*(-7*et**2-2+
        12*et*cosh(U)+(et**2-4)*cosh(2*U))*sinh(U)/(-1+et*cosh(U))**3)))

    if order<=1:
        return a0
    if order==2:
        return a0+a2
    if order==3:
        return a0+a2+a3


# ls=np.linspace(-l0,l0,1000)

# us0=get_u(ls,e1,eta,B1,0)
# us2=get_u(ls,e1,eta,B1,2)
# us3=get_u(ls,e1,eta,B1,3)



# fig, axs = plt.subplots(2, 1, sharex=True)
# # Remove horizontal space between axes
# fig.subplots_adjust(hspace=0.1)

# # Plot each graph, and manually set the y tick values
# axs[0].plot(ls,us0,label='Newtonian')
# axs[0].plot(ls,us2,label='2PN')
# axs[0].plot(ls,us3,label='3PN')

# axs[0].axhline(y=0,color='k', linestyle='--')
# axs[0].set_ylabel(r'$u$')
# axs[0].set_xlim([-l0,l0])


# axs[0].legend()



# axs[1].plot(ls,us2-us0,label=r'$\Delta_2 u$')
# axs[1].plot(ls,us3-us0,label=r'$\Delta_3 u$')

# axs[1].axhline(y=0,color='k', linestyle='--')

# axs[1].set_ylabel(r'$\Delta$')
# axs[1].legend()
# axs[1].set_xlabel(r'$l$')
# axs[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
# fig.suptitle('Mikkola Method:'"$\ e_{t}$= "+str(e1)+', b = '+str(B1)+r"$\ \frac{GM}{c^2}$")

# plt.savefig('mik.pdf')
# plt.show()
