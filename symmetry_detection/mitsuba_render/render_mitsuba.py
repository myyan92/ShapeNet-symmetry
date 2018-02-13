import mitsuba
import multiprocessing
import sys 
import os
import random
import numpy as np
from math import *
import shutil
import traceback

def normalize(modelName, newModelName, diameter = 1.0, YZswapping = 0):
    with open(modelName) as f:
        model = f.readlines()

    v = list()
    f = list()
    for line in model:
        if line.startswith("v "):
            line = line.strip().split()
            v.append([float(line[1]), float(line[2]), float(line[3]), 1])
        elif line.startswith("f "):
            line = line.strip().split()
            line = [l.split('/')[0] for l in line]
            f.append([int(line[1])-1, int(line[2])-1, int(line[3])-1])
    v = np.matrix(v).transpose()
    Tswap = np.matrix([[1.0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    T = np.matrix([[1.0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    if len(v) == 0 and YZswapping:
        return Tswap.tolist()
    if len(v) == 0 and YZswapping==0:
        return T.tolist()
    mv = np.zeros((3,1))
    totalarea = 0
    for face in f:
        area = np.linalg.norm(np.cross(v[0:3,face[0]]-v[0:3,face[1]], v[0:3,face[1]]-v[0:3,face[2]], axisa=0, axisb=0))
        mv += (v[0:3,face[0]] + v[0:3,face[1]] + v[0:3,face[2]]) * area / 3
        totalarea += area
    mv = mv / totalarea
    mx = mv[0, 0]
    my = mv[1, 0]
    mz = mv[2, 0]
    nv = v[0:3,:] - np.matrix([[mx],[my],[mz]])
    radius = float(np.max(np.sqrt((np.multiply(nv,nv)).sum(0))))
    T[0, 0] = diameter / radius
    T[1, 1] = diameter / radius
    T[2, 2] = diameter / radius
    T[0, 3] = -mx * diameter /radius
    T[1, 3] = -my * diameter /radius
    T[2, 3] = -mz * diameter /radius
    if YZswapping:
        T = Tswap.dot(T)
    v = T.dot(v)

    i = 0
    fout = open(newModelName, "w")
    for line in model:
        if line.startswith("v "):
            fout.write("v %f %f %f\n" % (v[0, i], v[1, i], v[2, i]))
            i += 1
        else:
            fout.write(line.strip() + "\n")
    fout.close()
    return T.tolist()


def render_model(model_name, synset):
    shutil.copytree(model_name, "./temp")

    md5 = model_name.split('/')[-1]
    print(md5)
    symmetry_file='../Results/' +synset+ '/'+md5+'.sym3t'
    if os.path.isfile(symmetry_file+'_g'):
        symmetry_file = symmetry_file+'_g'
    with open(symmetry_file) as f:
        lines = f.readlines()
    lines = [line.strip().split() for line in lines]
    sym_type = lines[0][0]
    if sym_type == "E":
        shutil.rmtree("./temp")
        return

#    translate = np.array([float(lines[1][0]), float(lines[1][1]), float(lines[1][2])])
#    right = np.array([float(lines[2][0]), float(lines[3][0]), float(lines[4][0])])
#    up = np.array([float(lines[2][1]), float(lines[3][1]), float(lines[4][1])])
#    front = np.array([float(lines[2][2]), float(lines[3][2]), float(lines[4][2])])

#    T = np.array([[1.0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
#    T[0, 0:3] = front / np.linalg.norm(front)
#    up = up - front * np.dot(front, up)
#    T[1, 0:3] = up / np.linalg.norm(up)
#    T[2, 0:3] = np.cross(up, front)
#    T[0:3, 3] = translate[:,np.newaxis]
#    print(T)

    normalize("./temp/model.obj", './temp/model.obj')
    transformParam = ' '.join([lines[2][2], lines[3][2], lines[4][2], '0', lines[2][1], lines[3][1], lines[4][1], '0', lines[2][0], lines[3][0], lines[4][0], '0', '0 0 0 1'])
    print(transformParam)
    os.system("sed 's,$mat,"+transformParam+",g' < render.xml > render_tmp.xml")
    os.system("sed -i 's,$x,"+lines[1][0]+",g' render_tmp.xml")
    os.system("sed -i 's,$y,"+lines[1][1]+",g' render_tmp.xml")
    os.system("sed -i 's,$z,"+lines[1][2]+",g' render_tmp.xml")
    scene = mitsuba.render.SceneHandler.loadScene(fileResolver.resolve("render_tmp.xml"))
    shape = scene.getShapes()
    print(shape)
    sensor = scene.getSensor()
    queue = mitsuba.render.RenderQueue()

    idx = 1
    if sym_type == "Cs":
        sensor_trans = [ mitsuba.core.Transform.lookAt(mitsuba.core.Point(6,2,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)),
                         mitsuba.core.Transform.lookAt(mitsuba.core.Point(6,-2,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)) ] 
        reflect_trans = mitsuba.core.Transform(mitsuba.core.Matrix4x4([1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0]))
        shape_trans = [mitsuba.core.Transform(), reflect_trans]
        for vi,sensor_t in enumerate(sensor_trans):
            for vj,shape_t in enumerate(shape_trans):
                sensor.setWorldTransform(mitsuba.core.Transform.inverse(shape_t)*sensor_t)
                output_path = './results/%s/%s_s%d_v%d.%d.png' % (synset, md5, idx, vi, vj) 
                scene.setDestinationFile(output_path)
                job = mitsuba.render.RenderJob('view_%d_%d_%d' % (idx, vi, vj), scene, queue)
                job.start()
                queue.waitLeft(0)
                queue.join()
    elif sym_type[0] in ['C', 'D']: #rotation
        try:
            degree = int(sym_type[1:])
        except:
            degree = int(sym_type[1:-1])
        sensor_trans = [ mitsuba.core.Transform.lookAt(mitsuba.core.Point(0,0,7),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)), 
                         mitsuba.core.Transform.lookAt(mitsuba.core.Point(0,7,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,0,1)) ] 
        shape_trans = [ mitsuba.core.Transform(),  mitsuba.core.Transform.rotate(mitsuba.core.Vector(0,1,0), float(360)/degree)]
        for vi,sensor_t in enumerate(sensor_trans):
            for vj,shape_t in enumerate(shape_trans):
                sensor.setWorldTransform(mitsuba.core.Transform.inverse(shape_t)*sensor_t)
                output_path = './results/%s/%s_s%d_v%d.%d.png' % (synset, md5, idx, vi, vj)
                scene.setDestinationFile(output_path)
                job = mitsuba.render.RenderJob('view_%d_%d_%d' % (idx, vi, vj), scene, queue)
                job.start()
                queue.waitLeft(0)
                queue.join()
        idx += 1
        if (sym_type[0]=='C' and sym_type[-1]=='v') or (sym_type[0]=='D' and sym_type[-1] in ['d','h']):
            for d in xrange(degree):
                sensor_trans = [ mitsuba.core.Transform.lookAt(mitsuba.core.Point(6,2,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)),
                                 mitsuba.core.Transform.lookAt(mitsuba.core.Point(6,-2,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)) ]
                reflect_trans = mitsuba.core.Transform(mitsuba.core.Matrix4x4([1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0]))
                rotation_trans =mitsuba.core.Transform.rotate(mitsuba.core.Vector(0,1,0), float(180)*d/degree) 
                shape_trans = [rotation_trans, reflect_trans*rotation_trans]
                for vi,sensor_t in enumerate(sensor_trans):
                    for vj,shape_t in enumerate(shape_trans):
                        sensor.setWorldTransform(mitsuba.core.Transform.inverse(shape_t)*sensor_t)
                        output_path = './results/%s/%s_s%d_v%d.%d.png' % (synset, md5, idx, vi, vj)
                        scene.setDestinationFile(output_path)
                        job = mitsuba.render.RenderJob('view_%d_%d_%d' % (idx, vi, vj), scene, queue)
                        job.start()
                        queue.waitLeft(0)
                        queue.join()
                idx += 1
        if sym_type[-1]=='h':
            sensor_trans = [ mitsuba.core.Transform.lookAt(mitsuba.core.Point(2,0,6),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)),
                             mitsuba.core.Transform.lookAt(mitsuba.core.Point(2,0,-6),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)) ]
            reflect_trans = mitsuba.core.Transform(mitsuba.core.Matrix4x4([1.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0]))
            shape_trans = [mitsuba.core.Transform(), reflect_trans]
            for vi,sensor_t in enumerate(sensor_trans):
                for vj,shape_t in enumerate(shape_trans):
                    sensor.setWorldTransform(mitsuba.core.Transform.inverse(shape_t)*sensor_t)
                    output_path = './results/%s/%s_s%d_v%d.%d.png' % (synset, md5, idx, vi, vj)
                    scene.setDestinationFile(output_path)
                    job = mitsuba.render.RenderJob('view_%d_%d_%d' % (idx, vi, vj), scene, queue)
                    job.start()
                    queue.waitLeft(0)
                    queue.join()
            idx += 1
        if sym_type[0]=='D':
            for d in xrange(degree):
                sensor_trans = [ mitsuba.core.Transform.lookAt(mitsuba.core.Point(7,0,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,1,0)),
                                 mitsuba.core.Transform.lookAt(mitsuba.core.Point(0,7,0),mitsuba.core.Point(0,0,0),mitsuba.core.Vector(0,0,1)) ]
                rotation_trans =mitsuba.core.Transform.rotate(mitsuba.core.Vector(0,1,0), float(180)*d/degree)
                shape_trans = [ rotation_trans,  mitsuba.core.Transform.rotate(mitsuba.core.Vector(1,0,0), float(180))*rotation_trans]
                for vi,sensor_t in enumerate(sensor_trans):
                    for vj,shape_t in enumerate(shape_trans):
                        sensor.setWorldTransform(mitsuba.core.Transform.inverse(shape_t)*sensor_t)
                        output_path = './results/%s/%s_s%d_v%d.%d.png' % (synset, md5, idx, vi, vj)
                        scene.setDestinationFile(output_path)
                        job = mitsuba.render.RenderJob('view_%d_%d_%d' % (idx, vi, vj), scene, queue)
                        job.start()
                        queue.waitLeft(0)
                        queue.join()
                idx += 1

    shutil.rmtree('./temp')

synset = sys.argv[1]
if not os.path.isdir('./results/' + synset):
    os.mkdir('./results/' + synset)
fileResolver = mitsuba.core.Thread.getThread().getFileResolver()
fileResolver.appendPath('./results/' + synset)

scheduler = mitsuba.core.Scheduler.getInstance()

#print(multiprocessing.cpu_count())
for i in range(0, 16): #multiprocessing.cpu_count()):
    scheduler.registerWorker(mitsuba.core.LocalWorker(i, 'wrk%i' % i))
scheduler.start()

filelist = '/orions4-zfs/projects/haohe/3DSIChallenge/deduplicate/deduplicated_model_list/'+synset+'.txt'
with open(filelist) as f:
    models = f.readlines()
models = [md5.strip() for md5 in models]
for i in range(len(models)):
    md5 = models[i]
#    md5="5334cfe15f80099f15ac67104577aee7"
    model_name = os.path.join('/orions3-zfs/projects/haosu/ShapeNetCore2015Spring/ShapeNetCore.v1', synset, 
                             md5)
    try:
        render_model(model_name, synset)
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        with open("errorLogs.txt", "a") as logfile:
            logfile.write("Processing %s" % (md5))
            logfile.write(''.join(lines))
        shutil.rmtree('./temp')
