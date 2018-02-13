import OpenEXR
import Imath
import cv2
import sys
import numpy as np

file = OpenEXR.InputFile(sys.argv[1])
pt = Imath.PixelType(Imath.PixelType.HALF)
window = file.header()['dataWindow']
size = (window.max.x - window.min.x + 1, window.max.y - window.min.y + 1)
#color pass
RGB = []
for c in ['B', 'G', 'R']:  #opencv is BGR
    str = file.channel('color.'+c, pt)
    img = np.fromstring(str, dtype=np.float16)
    img.shape = (size[1], size[0])
    img = np.expand_dims(img.astype(np.float32), axis=2)
    RGB.append(img)
RGB = np.concatenate(RGB, axis=2)
cv2.imshow('rgb', RGB)
cv2.waitKey(0)

#position pass
XYZ=[]
for c in ['B', 'G', 'R']:
    str = file.channel('position.'+c, pt)
    img = np.fromstring(str, dtype=np.float16)
    img.shape = (size[1], size[0])
    img = np.expand_dims(img.astype(np.float32), axis=2)
    XYZ.append(img)
XYZ = np.concatenate(XYZ, axis=2)
cv2.imshow('xyz', XYZ + 0.5)
cv2.waitKey(0)

