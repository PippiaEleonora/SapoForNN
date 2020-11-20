'''
@author: Nikos
'''

import numpy as np
import numpy.polynomial.polynomial as poly
#import torch
import matplotlib.pyplot as plt
from config import config

def poly_approx():
    upper =2
    lower=-2
    length=40
    x_poly=np.linspace(lower,upper,length)
    #x_all=[lower+x*(upper-lower)/length for x in range(length)]
    #x=torch.FloatTensor(x_all)
    #y=torch.tanh(x)

    y_poly=np.tanh(x_poly)
    coeffs=poly.polyfit(x_poly,y_poly,config.poly_order)
    for i in range(len(coeffs)):
        if abs(coeffs[i])<1e-7:
            coeffs[i]=0
    
    #print("\n The polynomial is a0+a1*x+a2*x^2+a3*x^3...\n\n",coeffs)
    #print("\n The polynomial order is {}.\n".format(config.poly_order))
    ffit=poly.polyval(x_poly,coeffs,".")
    plt.plot(x_poly,ffit,label="poly")
    plt.plot(x_poly,y_poly,marker="o",label="original",linestyle="none")
    #ffit=poly.Polynomial(coeffs)
    #plt.plot(x,ffit(x))
    plt.legend()
    plt.title("Testing polyfit against original")
    plt.xlabel("x")
    plt.ylabel("tanh(x)")
    plt.show(block=False)
    #print(type(coeffs))
    return coeffs
