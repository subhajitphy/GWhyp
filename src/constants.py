
import astropy.constants as ac

#Constants to be used later
GMsun = ac.GM_sun.value                                  #G*(Solar Mass)
c = ac.c.value                                           #Speed of light in vacuum
dsun = GMsun/(c**2)                                      #Natural length scale
tsun = GMsun/(c**3)                                      #Natural time scale
pc = ac.pc.value                                         #One parsec
yr = (365.25)*(24)*(60)*(60)        
