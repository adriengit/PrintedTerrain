# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 14:17:37 2015

3D Print of terrain elavation

@author: Adrien LE NAOUR
"""
import numpy as np
import matplotlib.pylab as plt
import os

class pystl():
    """ A Python class to manipulate 3D file format (.scad and .STL file format). """
    
    def read_hgt(self, filename):
        """ 
        Read .hgt file
        
        Examples
        ---------
        >>>plt.imshow(topo_read('datatopo/N48W005.hgt'))
        """
        siz = os.path.getsize(filename)
        dim = int(np.sqrt(siz/2))
        
        assert dim*dim*2 == siz, 'Invalid file size'
        
        return np.fromfile(filename, np.dtype('>i2'), dim*dim).reshape((dim, dim))


    def get_points_faces(self,s):
        """
        Points and Faces
        """    
        ntau = s.shape[1]
        nt = s.shape[0]
        fl = s.min() - 0.1 * np.abs(s.min() - s.max())
        nh = 2
        
        nptop = s.shape[0] * s.shape[1]
        nplow = nptop
        
        # POINTS
        npoints = nptop + nplow
        points = np.zeros((npoints, 3), dtype=int)
        
        for i in range(nt):
            for j in range(ntau):
                points[i * ntau + j, :] = np.array([i, j, s[i,j]])
        
        for i in range(nt):
            for j in range(ntau):
                points[nptop + i * ntau + j, :] = np.array([i, j, fl])
        
        # FACES
        nfaces = (nt-1)*(ntau-1)*2 + (nt-1)*(ntau-1)*2 + (nt-1)*(nh-1)*2 + (nt-1)*(nh-1)*2 + (ntau-1)*(nh-1)*2 + (ntau-1)*(nh-1)*2        
        faces = np.zeros((nfaces, 3))
        
        #top
        nfc = 0
        for j in range(nt-1):
            for i in range(ntau-1):
                faces[nfc, :] = [ntau * j  + i, ntau * j + i + 1, ntau * (j+1) + i]
                nfc +=1
                faces[nfc, :] = [ntau * j + i + 1,  ntau * (j+1) + i + 1, ntau * (j+1) + i]
                nfc += 1 
        #low
        for j in range(nt-1):
            for i in range(ntau-1):
                faces[nfc, :] = nptop + np.array([ntau * j  + i, ntau * (j+1) + i, ntau * j + i + 1])
                nfc +=1
                faces[nfc, :] = nptop + np.array([ntau * j + i + 1, ntau * (j+1) + i, ntau * (j+1) + i + 1])
                nfc +=1
        
        #front
        for i in range(ntau-1):
            faces[nfc, :] = [nptop + i, nptop + i + 1, i+1]
            nfc +=1
            faces[nfc, :] = [nptop + i, i + 1, i] 
            nfc +=1
        #back
        for i in range(ntau-1):
            faces[nfc, :] = ntau * (nt-1) + np.array([nptop + i, 1 + i, nptop + i + 1])
            nfc +=1
            faces[nfc, :] = ntau * (nt-1) + np.array([nptop + i, i, i + 1])
            nfc +=1
        
        #left
        for i in range(nt-1):
            faces[nfc, :] = np.array([ntau * i, ntau * (i + 1), nptop +  ntau * i])
            nfc +=1
            faces[nfc, :] = np.array([ntau * (i + 1), nptop +  ntau * (i+1), nptop +  ntau * i])
            nfc +=1
            
        #right
        for i in range(nt-1):
            faces[nfc, :] = np.array([ntau * i + ntau - 1, nptop +  ntau * i + ntau - 1, ntau * (i + 1) + ntau -1])
            nfc +=1
            faces[nfc, :] = np.array([nptop +  ntau * i + ntau - 1, nptop + ntau * (i+1) + ntau - 1, ntau * (i + 1) + ntau -1])
            nfc +=1
        return points, faces
    
    def _write_to_scad(points, faces, fn="poly.scad"):
        """
        write a polyhedron shape from points and faces list 
        """
        f = open(fn, "w+")
        f.write("polyhedron( points = [")
        for i in range(points.shape[0]):
            f.write("[ %f, %f, %f],"%(points[i, 0], points[i, 1], points[i, 2]))
        f.write("],")
        f.write("triangles= [")
        for i in range(faces.shape[0]):
            f.write("[ %d, %d, %d],"%(faces[i, 0], faces[i, 1], faces[i, 2]))
        f.write("]);")
        f.close()
    
    def surface_to_scad(self, s, filename):
        """
        Convert a surface to a scad polyhedron
        """
        points, faces = self.get_points_faces(s)        
        self._write_to_scad(points, faces, fn=filename)


    def _get_normal(self, pts):
        """
        Get normal vector from 3 tri-dimentional points pts.
        care to normal direction should be taken.
        """
        a = pts[1, :] - pts[0, :]
        b = pts[2, :] - pts[0, :]
        return np.cross(a, b)/np.sqrt(np.sum(np.abs(np.cross(a, b))**2))
    
    def _write_stl(self, points, faces, fn="menez.stl"):
        """
        surface to stl
        """
        f = open(fn, "w+")
        solidn = "menez"
        f.write("solid %s\n"%solidn)
        for i in range(faces.shape[0]):
            pt = points[faces[i].astype(int)]
            nor = self._get_normal(pt)
            f.write("facet normal %f %f %f\n"%(nor[0], nor[1], nor[2]))
            f.write("    outer loop\n")
            f.write("        vertex %e %e %e \n"%(pt[0, 0], pt[0, 1], pt[0, 2]))
            f.write("        vertex %e %e %e \n"%(pt[2, 0], pt[2, 1], pt[2, 2]))
            f.write("        vertex %e %e %e \n"%(pt[1, 0], pt[1, 1], pt[1, 2]))
            f.write("    endloop\n")
            f.write("endfacet\n")
        f.write("endsolid %s\n"%solidn)
        f.close()
    
    def surface_to_stl(self, s, filename):
        """
        Convert a surface to a stl
        """
        points, faces = self.get_points_faces(s)      
        self._write_stl(points, faces, fn=filename)

if __name__ == "__main__":
    mystl = pystl()
    #data = mystl.read_hgt("datatopo/N48W005.hgt")
    #plt.imshow(data)
    p = np.load('zone_aarberg.npy')
    #p = np.load('zone_thun.npy')
    s = np.arange(0, 30).reshape((10,3))
    s = np.ones((4,)).reshape((2,2))
    mystl.surface_to_stl(p[::5, ::5]/25., 'aarberg.stl')
    
    #mystl.surface_to_stl(data[600:1200, 600:800]/10., 'menez_quart.stl')
    
