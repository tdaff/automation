#!/usr/bin/env python

"""
gui_terminal

A frontend to look at structures in a character based terminal

"""

import curses
import time

import numpy as np
from numpy import dot

from faps import Structure


DEBUG = open('debug.txt', 'wb')
ASPECT_RATIO_Y_X = 1.6

def main(stdscr):
    """Wrapped function with all the programming stuff"""

    debug(curses.mousemask(curses.ALL_MOUSE_EVENTS))

    irmof1 = Structure("IMOF1")
    irmof1.from_cif("test.cif")#,
#                    cell=(25.8320, 25.832, 25.832, 90, 90, 90))

    zoom = 2*max(irmof1.cell.params[:3])
    origin = np.sum(irmof1.cell.cell, axis=0).tolist()

    origin = [0, 0, -50]
    w_vector = [0, 0, 1]
    u_vector = [1, 0, 0]
    v_vector = [0, 1, 0]
    direction = [0, 0, 1]
    up = [0, 1, 0]
    right = [1, 0, 0]

    orthographic(irmof1, stdscr, origin, zoom, direction, up, right)
    while 1:
        gotch = stdscr.getch()
        if gotch == ord('q'):
            break
        elif gotch == curses.KEY_LEFT:
            direction = dot(direction, rotation((0, 0, 10)))
            right = dot(right, rotation((0, 0, 10)))
            orthographic(irmof1, stdscr, origin, zoom, direction, up, right)
        elif gotch == curses.KEY_NPAGE:
            zoom += 1.0
            orthographic(irmof1, stdscr, origin, zoom, direction, up, right)
        elif gotch == curses.KEY_PPAGE:
            zoom -= 1.0
            orthographic(irmof1, stdscr, origin, zoom, direction, up, right)
        elif gotch == curses.KEY_MOUSE:
            debug(curses.getmouse())
        else:
            debug(curses.keyname(gotch))


def orthographic(struct, stdscr, origin, zoom, direction, up, right):
    """Put the structure on the screen."""
    stdscr.erase()
    maxy, maxx = stdscr.getmaxyx()
    debug(maxx, maxy)
    iatoms = [[atom.ipos(struct.cell.cell, struct.cell.inverse), atom] for atom in struct.atoms]
    atoms = [[[dot(right, atom[0]), dot(up, atom[0]), dot(direction, atom[0])], atom[1]] for atom in iatoms]
    for pos, atom in sorted(atoms, key=lambda x: x[0][2]):
        pos_x = int((maxx/2.0) - (maxx*((origin[1] - pos[0])/zoom)))
        pos_y = int(((maxx/2.0) - (maxx*((origin[0] - pos[1])/zoom)))/ASPECT_RATIO_Y_X)
        if 0 < pos_x < maxx and 0 < pos_y < maxy:
            stdscr.addstr(pos_y, pos_x, atom.type)
    stdscr.refresh()


def rotation(theta):
   tx,ty,tz = theta

   Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
   Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
   Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

   return np.dot(Rx, np.dot(Ry, Rz))

def rotation(theta, R = np.zeros((3,3))):
    cx,cy,cz = np.cos(theta)
    sx,sy,sz = np.sin(theta)
    R.flat = (cx*cz - sx*cy*sz, cx*sz + sx*cy*cz, sx*sy,
        -sx*cz - cx*cy*sz, -sx*sz + cx*cy*cz,
        cx*sy, sy*sz, -sy*cz, cy)
    return R

def make_axis_rotation_matrix(direction, angle):
     """
     Create a rotation matrix corresponding to the rotation around a general
     axis by a specified angle.

     R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)

     Parameters:

         angle : float a
         direction : array d
     """
     d = np.array(direction, dtype=np.float64)
     d /= np.linalg.norm(d)

     eye = np.eye(3, dtype=np.float64)
     ddt = np.outer(d, d)
     skew = np.array([[    0,  d[2],  -d[1]],
                      [-d[2],     0,  d[0]],
                      [d[1], -d[0],    0]], dtype=np.float64)

     mtx = ddt + np.cos(angle) * (eye - ddt) + np.sin(angle) * skew
     return mtx

def debug(*args):
    for arg in args:
        DEBUG.write(repr(arg) + '\n')
    DEBUG.flush()

curses.wrapper(main)
