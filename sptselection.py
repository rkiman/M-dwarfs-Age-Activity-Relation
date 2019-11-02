#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:33:55 2019

@author: rociokiman
"""

import numpy as np


spt_n1 = ['F0','F0.1','F0.2','F0.3','F0.4','F0.5','F0.6','F0.7','F0.8','F0.9',
          'F1','F1.1','F1.2','F1.3','F1.4','F1.5','F1.6','F1.7','F1.8','F1.9',
          'F2','F2.1','F2.2','F2.3','F2.4','F2.5','F2.6','F2.7','F2.8','F2.9',
          'F3','F3.1','F3.2','F3.3','F3.4','F3.5','F3.6','F3.7','F3.8','F3.9',
          'F4','F4.1','F4.2','F4.3','F4.4','F4.5','F4.6','F4.7','F4.8','F4.9',
          'F5','F5.1','F5.2','F5.3','F5.4','F5.5','F5.6','F5.7','F5.8','F5.9',
          'F6','F6.1','F6.2','F6.3','F6.4','F6.5','F6.6','F6.7','F6.8','F6.9',
          'F7','F7.1','F7.2','F7.3','F7.4','F7.5','F7.6','F7.7','F7.8','F7.9',
          'F8','F8.1','F8.2','F8.3','F8.4','F8.5','F8.6','F8.7','F8.8','F8.9',
          'F9','F9.1','F9.2','F9.3','F9.4','F9.5','F9.6','F9.7','F9.8','F9.9',
          'G0','G0.1','G0.2','G0.3','G0.4','G0.5','G0.6','G0.7','G0.8','G0.9',
          'G1','G1.1','G1.2','G1.3','G1.4','G1.5','G1.6','G1.7','G1.8','G1.9',
          'G2','G2.1','G2.2','G2.3','G2.4','G2.5','G2.6','G2.7','G2.8','G2.9',
          'G3','G3.1','G3.2','G3.3','G3.4','G3.5','G3.6','G3.7','G3.8','G3.9',
          'G4','G4.1','G4.2','G4.3','G4.4','G4.5','G4.6','G4.7','G4.8','G4.9',
          'G5','G5.1','G5.2','G5.3','G5.4','G5.5','G5.6','G5.7','G5.8','G5.9',
          'G6','G6.1','G6.2','G6.3','G6.4','G6.5','G6.6','G6.7','G6.8','G6.9',
          'G7','G7.1','G7.2','G7.3','G7.4','G7.5','G7.6','G7.7','G7.8','G7.9',
          'G8','G8.1','G8.2','G8.3','G8.4','G8.5','G8.6','G8.7','G8.8','G8.9',
          'G9','G9.1','G9.2','G9.3','G9.4','G9.5','G9.6','G9.7','G9.8','G9.9',
          'K0','K0.1','K0.2','K0.3','K0.4','K0.5','K0.6','K0.7','K0.8','K0.9',
          'K1','K1.1','K1.2','K1.3','K1.4','K1.5','K1.6','K1.7','K1.8','K1.9',
          'K2','K2.1','K2.2','K2.3','K2.4','K2.5','K2.6','K2.7','K2.8','K2.9',
          'K3','K3.1','K3.2','K3.3','K3.4','K3.5','K3.6','K3.7','K3.8','K3.9',
          'K4','K4.1','K4.2','K4.3','K4.4','K4.5','K4.6','K4.7','K4.8','K4.9',
          'K5','K5.1','K5.2','K5.3','K5.4','K5.5','K5.6','K5.7','K5.8','K5.9',
          'K6','K6.1','K6.2','K6.3','K6.4','K6.5','K6.6','K6.7','K6.8','K6.9',
          'K7','K7.1','K7.2','K7.3','K7.4','K7.5','K7.6','K7.7','K7.8','K7.9',
          'K8','K8.1','K8.2','K8.3','K8.4','K8.5','K8.6','K8.7','K8.8','K8.9',
          'K9','K9.1','K9.2','K9.3','K9.4','K9.5','K9.6','K9.7','K9.8','K9.9',
          'M0','M0.1','M0.2','M0.3','M0.4','M0.5','M0.6','M0.7','M0.8','M0.9',
          'M1','M1.1','M1.2','M1.3','M1.4','M1.5','M1.6','M1.7','M1.8','M1.9',
          'M2','M2.1','M2.2','M2.3','M2.4','M2.5','M2.6','M2.7','M2.8','M2.9',
          'M3','M3.1','M3.2','M3.3','M3.4','M3.5','M3.6','M3.7','M3.8','M3.9',
          'M4','M4.1','M4.2','M4.3','M4.4','M4.5','M4.6','M4.7','M4.8','M4.9',
          'M5','M5.1','M5.2','M5.3','M5.4','M5.5','M5.6','M5.7','M5.8','M5.9',
          'M6','M6.1','M6.2','M6.3','M6.4','M6.5','M6.6','M6.7','M6.8','M6.9',
          'M7','M7.1','M7.2','M7.3','M7.4','M7.5','M7.6','M7.7','M7.8','M7.9',
          'M8','M8.1','M8.2','M8.3','M8.4','M8.5','M8.6','M8.7','M8.8','M8.9',
          'M9','M9.1','M9.2','M9.3','M9.4','M9.5','M9.6','M9.7','M9.8','M9.9']

spt_list_num1 = np.arange(-300,100)/10


spt_n2 = ['O0','O1','O2','O3','O4','O5','O6','O7','O8','O9',
          'B0','B1','B2','B3','B4','B5','B6','B7','B8','B9',
          'A0','A1','A2','A3','A4','A5','A6','A7','A8','A9',
          'F0','F1','F2','F3','F4','F5','F6','F7','F8','F9',
          'G0','G1','G2','G3','G4','G5','G6','G7','G8','G9',
          'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9',
          'M0','M1','M2','M3','M4','M5','M6','M7','M8','M9',
          'L0','L1','L2','L3','L4','L5','L6','L7','L8','L9']

spt_list_num2 = np.arange(-60,20)

spt_n3 = ['K5V','K5.5V','K6V','K6.5V','K7V','K7.5V','K8V','K8.5V','K9V','K9.5V',
          'M0V','M0.5V','M1V','M1.5V','M2V','M2.5V','M3V','M3.5V','M4V','M4.5V',
          'M5V','M5.5V','M6V','M6.5V','M7V','M7.5V','M8V','M8.5V','M9V','M9.5V']

spt_list_num3 = np.arange(-10,20)/2

spt_n4 = ['K5Ve','K5.5Ve','K6Ve','K6.5Ve','K7Ve','K7.5Ve','K8Ve','K8.5Ve',
          'K9Ve','K9.5Ve',
          'M0Ve','M0.5Ve','M1Ve','M1.5Ve','M2Ve','M2.5Ve','M3Ve','M3.5Ve',
          'M4Ve','M4.5Ve','M5Ve','M5.5Ve','M6Ve','M6.5Ve','M7Ve','M7.5Ve',
          'M8Ve','M8.5Ve','M9Ve','M9.5Ve']

spt_list_num4 = np.arange(-10,20)/2

spt_n42 = ['M0.0Ve','M1.0Ve','M2.0Ve','M3.0Ve','M4.0Ve','M5.0Ve','M6.0Ve',
           'M7.0Ve','M8.0Ve','M9.0Ve']

spt_list_num42 = np.arange(0,10)

spt_n43 = ['K5V ','K5.5V ','K6V ','K6.5V ','K7V ','K7.5V ','K8V ','K8.5V ',
          'K9V ','K9.5V ',
          'M0V ','M0.5V ','M1V ','M1.5V ','M2V ','M2.5V ','M3V ','M3.5V ',
          'M4V ','M4.5V ','M5V ','M5.5V ','M6V ','M6.5V ','M7V ','M7.5V ',
          'M8V ','M8.5V ','M9V ','M9.5V ']

spt_list_num43 = np.arange(-10,20)/2

spt_n44 = ['M0.0V ','M1.0V ','M2.0V ','M3.0V ','M4.0V ','M5.0V ','M6.0V ',
           'M7.0V ','M8.0V ','M9.0V ']

spt_list_num44 = np.arange(0,10)

spt_n5 = ['G0.0','G1.0','G2.0','G3.0','G4.0','G5.0','G6.0','G7.0','G8.0','G9.0',
          'K0.0','K1.0','K2.0','K3.0','K4.0','K5.0','K6.0','K7.0','K8.0','K9.0',
          'M0.0','M1.0','M2.0','M3.0','M4.0','M5.0','M6.0','M7.0','M8.0','M9.0']

spt_list_num5 = np.arange(-20,10)

spt_list_name = np.array(spt_n1 + spt_n2 + spt_n3 + spt_n4 + spt_n5 + spt_n42 + 
                         spt_n43 + spt_n44) 
spt_list_num = np.concatenate((spt_list_num1,spt_list_num2,
                               spt_list_num3,spt_list_num4,
                               spt_list_num5,spt_list_num42,
                               spt_list_num43,spt_list_num44))

def organize_spt(spt_sample):
    spt_sample_good = []
    for x in spt_sample:
        #print(x)
        if(str(x)=='G/K '):
            spt_sample_good.append(np.nan)
        elif(str(x)=='M0   2'):
            spt_sample_good.append(0)
        elif(str(x)=='M8 J  '):
            spt_sample_good.append(8)
        elif(str(x)=='M6.25'):
            spt_sample_good.append(6.25)
        elif(str(x)=='M4.75'):
            spt_sample_good.append(4.75)
        elif(str(x)==''):
            spt_sample_good.append(np.nan)
        elif(str(x)=='G10'):
            spt_sample_good.append(np.nan)
        elif(str(x)=='F10'):
            spt_sample_good.append(np.nan)
        elif(str(x)=='A10'):
            spt_sample_good.append(np.nan)
        elif(str(x)=='n/a'):
            spt_sample_good.append(np.nan)
        elif(str(x)=='M0V(e)'):
            spt_sample_good.append(0.0)
        elif(str(x)=='M0IIIe'):
            spt_sample_good.append(0.0)
        elif(str(x)=='M3V(e)'):
            spt_sample_good.append(3.0)
        elif(str(x)=='M2V(e)'):
            spt_sample_good.append(2.0)
        elif(str(x)=='M2III'):
            spt_sample_good.append(2.0)
        elif(str(x)=='M4IVe'):
            spt_sample_good.append(4.0)
        elif(str(x)=='M5IVe'):
            spt_sample_good.append(5.0)
        elif(str(x)=='M5III'):
            spt_sample_good.append(5.0)
        elif(str(x)=='M0e'):
            spt_sample_good.append(0.0)
        elif('+' in list(x)):
            spt_sample_good.append(np.nan)
        elif('-' in list(x)):
            spt_sample_good.append(np.nan)
        elif('*' in list(x)):
            spt_sample_good.append(np.nan)
        elif(('s' in list(x)) and ('b' in list(x))):
            spt_sample_good.append(np.nan)
        elif(('v' in list(x)) and ('b' in list(x))):
            spt_sample_good.append(np.nan)
        else:
            mask = str(x) == spt_list_name
            spt_sample_good.append(spt_list_num[mask][0])
    return np.array(spt_sample_good) 