#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits 


# In[3]:


def load_file(files):
    if type(files) != str:
        print('input needs to be a string')
        return
    if (files.endswith('.fits')) | (files.endswith('.Fits') | files.endswith('.fit') | files.endswith('.Fit')):
        hdu = fits.open(files)
        header = hdu[0].header
        data = hdu[0].data
        return header, data
    else:
        print('incorrect file extension: please use [.fits] [.Fits] [.fit] [.Fit] files')


# In[4]:


def implot(dat, **kwargs):
    if ('fig') in kwargs:
        fig = kwargs['fig']
    fig = plt.figure(figsize = (15, 13))

    cmap = 'gray_r'
    mu = np.mean(dat)
    sigma = np.std(dat)
    scale = .5
    if ('scale') in kwargs:
        scale = kwargs['scale']
    alpha = mu - scale*sigma
    beta = mu + scale*sigma
    vmin = alpha
    vmax = beta
    if ('vmin') in kwargs:
        vmin = kwargs['vmin']
    if ('vmax') in kwargs:
        vmax = kwargs['vmax']
    if ('cmap') in kwargs:
        cmap = kwargs['cmap']

    plt.imshow(dat, vmin = vmin, vmax = vmax, origin = 'lower', cmap = cmap)
    
    plt.show()


# In[ ]:




