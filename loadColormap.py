# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 23:07:20 2013

@author: Reinhardt A.W. Maier <rema@zaehlwerk.net>
"""
def whiteCMap():
    from numpy import genfromtxt
    from matplotlib.colors import ListedColormap
    
    whiteColormap = genfromtxt('colormap.txt', delimiter=',')
    whiteCMap = ListedColormap(whiteColormap, name='white')
    
    return whiteCMap
