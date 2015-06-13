from numpy import *
from pylab import *

data = loadtxt('5p1s6e/fields.data', type('string'))

l = transpose(data)[1]
b = transpose(data)[2]
R = transpose(data)[3]
l = map(float, l)
b = map(float, b)
R = map(float, R)

def xy(l, b, w = 4096):
    n = -1.0 if (b < 0) else +1.
    x = w/2.0 * sqrt(1 - n*sin(radians(b))) * cos(radians(l)) + (w/2.0 - 0.5)
    y = -w/2.0 * n * sqrt(1 - n*sin(radians(b))) * sin(radians(l)) + (w/2.0 - 0.5)
    return [x,y]
    
w = 4096.0



xs = []
ys = []
rs = []

xn = []
yn = []
rn = []

datas = zeros((w, w))
datan = zeros((853, 853))
datan = datas

for i in range(0, len(R)):
    bi = b[i]
    li = l[i]
    Ri = R[i]
    if (bi < 0):
        # SGP
        n = -1.0
        x = w/2.0 * sqrt(1 - n*sin(radians(bi))) * cos(radians(li)) + (w/2.0 - 0.5)
        y = -w/2.0 * n * sqrt(1 - n*sin(radians(bi))) * sin(radians(li)) + (w/2.0 - 0.5)
        datas[len(xs)][len(ys)] = Ri
        xs.append(x)
        ys.append(y)
        rs.append(Ri)
    
    else:
        # NGP
        n = 1.0
        x = w/2.0 * sqrt(1 - n*sin(radians(bi))) * cos(radians(li)) + (w/2.0 - 0.5)
        y = -w/2.0 * n * sqrt(1 - n*sin(radians(bi))) * sin(radians(li)) + (w/2.0 - 0.5)

        #datan[len(xn)][len(yn)] = Ri
        datan[round(x)][round(y)] = Ri
        xn.append(x)
        yn.append(y)
        rn.append(Ri)



subplot(121)
#X, Y = meshgrid(xn, yn)

imshow(datan, interpolation='bilinear', cmap=gray, origin='lower')
drange = range(0, int(w))
#contour(xn, yn, datan, interpolaton='bilinear')
colorbar()
#Test projection
bi = 30.0
li = 30.0
l_test = arange(0, 360, li)
b_test = arange(0, 90, bi)

for l_t in l_test:
    x, y = xy(l_t, range(0, 90))
    c = tuple([l_t/360.]*3)
    c = (0.5, 0.5, 0.5)
    plot(x,y, '--', c=c)
    
for b_t in b_test:
    x, y = xy(range(0, 360), b_t)
    c = tuple([b_t/90.]*3)
    c = (0.5, 0.5, 0.5)
    plot(x,y, '--', c=c)

subplot(122)
#contour(xs, ys, rs)
#contour(xs, ys, datas)
#Test projection
bi = 30.0
li = 30.0
l_test = arange(0, 360, li)
b_test = arange(0, 90, bi)

for l_t in l_test:
    x, y = xy(l_t, range(0, 90))
    c = tuple([l_t/360.]*3)
    c = (0.5, 0.5, 0.5)
    plot(x,y, '--', c=c)
    
    
for b_t in b_test:
    x, y = xy(range(0, 360), b_t)
    c = tuple([b_t/90.]*3)
    c = (0.5, 0.5, 0.5)
    plot(x,y, '--', c=c)

show()