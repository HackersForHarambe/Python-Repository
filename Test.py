import pygame
import cmath
import math as m
from pygame import *
from cmath import *
import matplotlib as mat
import numpy as np


xmax = 20
ymax = 16
inter = 1
interx = 1
tmax = 10
intert = 10

j = sqrt(-1)

WOUT_Last = [ 0, 0]


def ln(x):
    try:
        return log(abs(x)) - j * ( atan( x.real / ( .0000001 + x.imag)) - pi / 2)
    except:
        return -100

def Integrate(f, x, D):
    Intout = 0
    for i in range( -int( D * abs(x)), int( D * abs(x))):
        Intout += f(i / D * (x + .00001) / (abs(x) + .00001)) / D
    return Intout

def Dedekind(x):
    DedeOUT = x
    if abs(x) <= 1:
        for i in range(1,101):
            DedeOUT *= ( 1 - x ** i)
    else:
        DedeOUT = 0
    return x * DedeOUT

def Zeta(z):
    ZetaOUT = 0
    if z.real == 1 and z.imag == 0:
        return 10000
    elif z.real >= 0:
        for i in range(1,501):
            ZetaOUT += (-1) ** i / (i ** z)
        return ZetaOUT / (2 ** (1 - z) - 1)
    else:
        return 0

def Eta(x):
    return Zeta(j * x + 1 / 2)

def RiemannZ(n):
    i = 0
    if n > 0:
        while int(n) != 0:
            while abs(Eta(i)) > .04:
                i += abs(Eta(i)) ** 2 / 5
            n += -1
            i += .1
    else:
        while int(n) != 0:
            while abs(Eta(-i)) > .04:
                i += -abs(Eta(i)) ** 2 / 5
            n += 1
            i += -.1
    return i

def fact(x):
    factOUT = 1
    for i in range(1, x + 1):
        factOUT *= i
    if factOUT == 1:
        return 1
    else:
        return factOUT

def liPart(k):
    liPOUT = 0
    for n in range(0, m.floor((k + 1) / 2)):
        liPOUT += 1/( 2 * n + 1)
    return liPOUT

def Welliptic(x, w1, w2):
    global WOUT2
    global WOUT_Last

    WOUT1 = 0
    for m_1 in range(-50,51):
        try:
            WOUT1 += sin( pi * ( x - 2 * w2 * m_1) / ( 2 * w1)) ** (-2)
        except:
            WOUT1 += 100
    if WOUT_Last != [ w1, w2]:
        WOUT2 = 0
        for m_2 in range(-50, 51):
            if m_2 == 0:
                continue
            else:
                WOUT2 += sin(pi * w2 * m_2 / w1) ** (-2)
    WOUT_Last = [w1, w2]
    return (pi / (2 * w1)) ** 2 * ( -1 / 3 + WOUT1 - WOUT2)

def li( z, p = 1):
    liOUT = 0
    for i in range(1, 100):
        liOUT += ( p ** i) * ( ln(z) ** i) * liPart(i) / ( fact(i) * (-2) ** ( i -1))
    if z.real == 1 and z.imag == 0:
        return -100
    else:
        return sqrt(z ** p) * liOUT + ln(abs(ln(z))) + ln(p) + np.euler_gamma

def Draw(F):
    for i_1 in range( -400, 401, interx):
        i_x = i_1 * xmax / 800
        draw.line(DisplaySurf, (0, 255, 0), (i_1 + 400, 300 - F(i_x).real * 600 / ymax),
            (i_1 + interx + 400, 300 - F(i_x + interx * xmax / 800).real * 600 / ymax), 3)
        draw.line(DisplaySurf, (255, 0, 0), (i_1  + 400, 300 - F(i_x).imag * 600 / ymax),
            (i_1 + interx + 400, 300 - F(i_x + interx * xmax / 800).imag * 600 / ymax), 3)

def DrawP(G, F):
    for i_1 in range( -tmax * intert, tmax * intert, 1):
        i_t = i_1 / intert
        draw.line(DisplaySurf, (0, 255, 0), ( 400 + G(i_t).real * 400 / xmax, 300 + F(i_t).real * 300 / ymax),
            ( 400 + G(i_t + 1 / intert).real * 400 / xmax,
              300 + F(i_t + 1 / intert).real * 300 / ymax), 3)
        draw.line(DisplaySurf, (255, 0, 0), ( 400 + G(i_t).imag * 400 / xmax, 300 + F(i_t).real * 300 / ymax),
            ( 400 + G(i_t + 1 / intert).imag * 400 / xmax,
              300 + F(i_t + 1 / intert).real * 300 / ymax), 3)

def Cloth(x):
    if x == 0:
        return 0
    else:
        return tan(x ** 2) + Cloth( x - 1)

init()
DisplaySurf = pygame.display.set_mode((800, 600), 0, 32)
display.set_caption('Graph')

DisplaySurf.fill(( 255, 255, 255))


for i_1 in range( -2 * (xmax), 2 * (xmax + inter), inter):
    draw.line(DisplaySurf, (51, 51, 51), (i_1 / 2 * 400 / xmax + 400, 600), (i_1 / 2 * 400 / xmax + 400, 0), 1)

for i_1 in range( -xmax - xmax % inter, xmax + inter, 2 * inter):
    draw.line(DisplaySurf, (102, 102, 102), ( i_1 * 400 / xmax + 400, 600), ( i_1 * 400 / xmax + 400, 0), 1)

for i_1 in range( -2 * (ymax), 2 * (ymax + inter), inter):
    draw.line(DisplaySurf, (51, 51, 51), ( 0, i_1 / 2 * 300 / ymax + 300), ( 800, i_1 / 2 * 300 / ymax + 300), 1)

for i_1 in range( -ymax - ymax % 2, ymax + inter, 2 * inter):
    draw.line(DisplaySurf, (102, 102, 102), ( 0, i_1 * 300 / ymax + 300), ( 800, i_1 * 300 / ymax + 300), 1)

draw.line(DisplaySurf, (0, 0, 0), ( 0, 300), ( 800, 300), 3)
draw.line(DisplaySurf, (0, 0, 0), ( 400, 0), ( 400, 600), 3)
draw.circle(DisplaySurf, ( 0, 0, 0), ( 400, 300), 8, 2)


Draw(lambda x: Integrate( lambda t: t ** 2 * cos(2 * pi * t * int(x)), 30, 10) + j * 5)


while True:
    for event in pygame.event.get():
        if event.type == QUIT:
            quit()
            exit()
    display.update()